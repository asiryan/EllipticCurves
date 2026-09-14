using System.Diagnostics;
using System.Text.Json;

// Parameter selection only. Finite-field scores never certify a rank.
sealed record SampleOptions(int Samples, int Height, int Keep, int RefineKeep,
    int FinalKeep, int PrimeBound, int Workers, int Seed);
sealed class SampleCheckpoint
{
    public SampleOptions Options { get; set; } = null!;
    public ulong RandomState { get; set; }
    public long Draws { get; set; }
    public long PrimitiveDraws { get; set; }
    public Candidate[] Retained { get; set; } = [];
}

static class CandidateSampler
{
    const int Scale = 1 << 20;
    sealed class Table
    {
        public readonly int P;
        readonly int[] inverse, scores;
        public Table(int p)
        {
            P = p; inverse = new int[p]; scores = new int[p + 1]; inverse[1] = 1;
            for (int i = 2; i < p; i++) inverse[i] = p - (int)((long)(p / i) * inverse[p % i] % p);
            var local = new Local302(p);
            for (int i = 0; i <= p; i++) scores[i] = (int)Math.Round(Scale * local.Score(i), MidpointRounding.AwayFromZero);
        }
        public int Score(int u, int v)
        {
            int d = v % P;
            return scores[d == 0 ? P : (int)((long)Local302.Mod(u, P) * inverse[d] % P)];
        }
    }
    static ulong Next(ref ulong state)
    {
        ulong z = unchecked(state += 0x9E3779B97F4A7C15UL);
        z = unchecked((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL);
        z = unchecked((z ^ (z >> 27)) * 0x94D049BB133111EBUL);
        return z ^ (z >> 31);
    }
    static int Gcd(int a, int b) { a = Math.Abs(a); while (b != 0) (a, b) = (b, a % b); return a; }
    static Candidate[] Ordered(IEnumerable<Candidate> rows, int count) => rows
        .OrderByDescending(c => c.Score).ThenBy(c => c.V).ThenBy(c => c.U).Take(count).ToArray();

    sealed class BandCheckpoint
    {
        public SampleOptions Options { get; set; } = null!;
        public ulong RandomState { get; set; }
        public long Draws { get; set; }
        public long PrimitiveDraws { get; set; }
        public Candidate[][] Retained { get; set; } = [[], [], [], []];
    }

    public static void RunBands(SampleOptions options, string directory, CancellationToken token)
    {
        Directory.CreateDirectory(directory);
        string checkpoint = Path.Combine(directory, "band-checkpoint.json");
        var state = File.Exists(checkpoint)
            ? JsonSerializer.Deserialize<BandCheckpoint>(File.ReadAllText(checkpoint))!
            : new BandCheckpoint { Options = options, RandomState = (ulong)options.Seed };
        if (state.Options != options) throw new ArgumentException("Band sampling options differ from checkpoint.");
        if (File.Exists(Path.Combine(directory, "complete.json"))) return;
        var watch = Stopwatch.StartNew();
        int[] Spread(int low, int high, int count)
        {
            var available = Local302.Primes(high).Where(p => p > low).ToArray();
            return Enumerable.Range(0, count).Select(i => available[(int)((long)(2*i+1)*available.Length/(2*count))]).Distinct().ToArray();
        }
        // Independent channels: a poor small-prime score cannot veto the later bands.
        // The fourth channel is a fixed-size random reservoir, independent of scores.
        var bands = new[] { Local302.Primes(1021), Spread(1021,8191,12), Spread(8191,32749,8) };
        var tables = new Table[bands.Length][];
        for (int band = 0; band < bands.Length; band++)
        {
            tables[band] = new Table[bands[band].Length]; int current = band;
            Parallel.For(0, bands[band].Length,
                new ParallelOptions { MaxDegreeOfParallelism = options.Workers, CancellationToken = token },
                i => tables[current][i] = new Table(bands[current][i]));
            Console.WriteLine($"band tables {band+1}/{bands.Length}: {watch.Elapsed.TotalSeconds:F2}s");
        }
        int quota = Math.Max(1, options.Keep/4);
        var queues = Enumerable.Range(0,4).Select(_ => new PriorityQueue<Candidate,(double,int,int)>()).ToArray();
        var members = Enumerable.Range(0,4).Select(_ => new HashSet<(int,int)>()).ToArray();
        void Retain(int channel, Candidate row)
        {
            if (members[channel].Contains((row.U,row.V))) return;
            var priority = (row.Score,-row.V,-row.U);
            if (queues[channel].Count >= quota)
            {
                queues[channel].TryPeek(out _,out var worst);
                if (priority.CompareTo(worst) <= 0) return;
                var removed = queues[channel].Dequeue(); members[channel].Remove((removed.U,removed.V));
            }
            queues[channel].Enqueue(row,priority); members[channel].Add((row.U,row.V));
        }
        for (int channel=0;channel<4;channel++) foreach(var row in state.Retained[channel]) Retain(channel,row);
        ulong randomState = state.RandomState;
        void Checkpoint()
        {
            state.RandomState=randomState;
            state.Retained=queues.Select(q=>Ordered(q.UnorderedItems.Select(e=>e.Element),quota)).ToArray();
            Hunt.Save(checkpoint,state);
        }
        try
        {
            while(state.PrimitiveDraws<options.Samples && state.Draws<8L*options.Samples)
            {
                token.ThrowIfCancellationRequested();state.Draws++;
                int u=(int)(Next(ref randomState)%(uint)(2L*options.Height+1))-options.Height;
                int v=1+(int)(Next(ref randomState)%(uint)options.Height);
                if(Gcd(u,v)!=1)continue;
                state.PrimitiveDraws++;
                double randomPriority=(Next(ref randomState)>>11)*(1.0/(1UL<<53));
                if(u!=164518 || v!=924945)
                {
                    for(int channel=0;channel<3;channel++)
                    {
                        long score=0;foreach(var table in tables[channel])score+=table.Score(u,v);
                        Retain(channel,new Candidate{U=u,V=v,Score=(double)score/Scale});
                    }
                    Retain(3,new Candidate{U=u,V=v,Score=randomPriority});
                }
                if(state.PrimitiveDraws%100000==0)
                {Checkpoint();Console.WriteLine($"band sample {state.PrimitiveDraws}/{options.Samples}: {watch.Elapsed.TotalSeconds:F2}s");}
            }
        }
        finally {Checkpoint();}
        var candidates=state.Retained.SelectMany(x=>x).DistinctBy(c=>(c.U,c.V))
            .Select(c=>new Candidate{U=c.U,V=c.V}).ToList();
        // All channel survivors see the full later interval: no intervening
        // cumulative-score cutoff (which rejected the rank-31 control).
        Hunt.ScoreStage(candidates,0,16381,options.Workers,token);
        foreach(var c in candidates)c.ScreeningScore=c.Score;
        Hunt.ScoreStage(candidates,16381,options.PrimeBound,options.Workers,token);
        foreach(var c in candidates)c.ValidationScore=c.Score-c.ScreeningScore;
        var full=Ordered(candidates,candidates.Count);
        var tail=candidates.OrderByDescending(c=>c.ValidationScore).ThenBy(c=>c.V).ThenBy(c=>c.U).ToArray();
        var final=new List<Candidate>();var seen=new HashSet<(int,int)>();
        for(int i=0;i<candidates.Count && final.Count<options.FinalKeep;i++)
            foreach(var c in new[]{full[i],tail[i]})
                if(final.Count<options.FinalKeep && seen.Add((c.U,c.V)))final.Add(c);
        var seenJ=new HashSet<string>{Family302.Curve(164518,924945).JInvariant.ToString()};
        var output=new List<object>();
        foreach(var c in final)
        {
            token.ThrowIfCancellationRequested();var (curve,points)=Structured302.Create(c.U,c.V);
            string j=curve.JInvariant.ToString();if(!seenJ.Add(j))continue;
            var proof=curve.GetRankLowerBound(points,1009,token);string file=$"curve_{c.U}_{c.V}.json";
            Hunt.Save(Path.Combine(directory,file),Hunt.CurveData(curve,points,$"icarm302-17 sections at {c.U}/{c.V}",proof.LowerBound));
            output.Add(new{id=$"{c.U}_{c.V}",u=c.U,v=c.V,file,score=c.Score,screening_score=c.ScreeningScore,
                tail_score=c.ValidationScore,seed_lower_bound=proof.LowerBound,j_invariant=j,published_novelty_verified=false});
        }
        Hunt.Save(Path.Combine(directory,"candidates.json"),output);
        Hunt.Save(Path.Combine(directory,"complete.json"),new{options,selection="independent bands plus random reservoir",
            bands,state.Draws,state.PrimitiveDraws,rescored=candidates.Count,exported=output.Count,
            seconds=watch.Elapsed.TotalSeconds,score_is_not_a_rank_bound=true,
            full_rescore_uses_screening_primes=true,unbiased_holdout_claimed=false});
        Console.WriteLine($"Prepared {output.Count} band-selected candidates in {watch.Elapsed.TotalSeconds:F2}s.");
    }

    public static void Run(SampleOptions options, string directory, CancellationToken token)
    {
        Directory.CreateDirectory(directory);
        string checkpointPath = Path.Combine(directory, "sample-checkpoint.json");
        var state = File.Exists(checkpointPath)
            ? JsonSerializer.Deserialize<SampleCheckpoint>(File.ReadAllText(checkpointPath))!
            : new SampleCheckpoint { Options = options, RandomState = (ulong)options.Seed };
        if (state.Options != options) throw new ArgumentException("Sampling options differ from checkpoint.");
        if (File.Exists(Path.Combine(directory, "complete.json"))) { Console.WriteLine("Candidate pool already complete."); return; }
        var timer = Stopwatch.StartNew();
        var queue = new PriorityQueue<Candidate, (double, int, int)>();
        var retained = new HashSet<(int, int)>();
        void Push(Candidate c) { queue.Enqueue(c, (c.Score, -c.V, -c.U)); retained.Add((c.U, c.V)); }
        foreach (var c in state.Retained) Push(c);
        var primes = Local302.Primes(1021); var tables = new Table[primes.Length];
        Parallel.For(0, primes.Length, new ParallelOptions { MaxDegreeOfParallelism = options.Workers, CancellationToken = token },
            i => tables[i] = new Table(primes[i]));
        ulong randomState = state.RandomState;
        void Checkpoint()
        {
            state.RandomState = randomState;
            state.Retained = Ordered(queue.UnorderedItems.Select(x => x.Element), options.Keep);
            Hunt.Save(checkpointPath, state);
        }
        try
        {
            while (state.PrimitiveDraws < options.Samples && state.Draws < 8L * options.Samples)
            {
                token.ThrowIfCancellationRequested(); state.Draws++;
                int u = (int)(Next(ref randomState) % (uint)(2L * options.Height + 1)) - options.Height;
                int v = 1 + (int)(Next(ref randomState) % (uint)options.Height);
                if (Gcd(u, v) != 1) continue;
                state.PrimitiveDraws++;
                // The known record is excluded, never injected as a candidate.
                if ((u != 164518 || v != 924945) && !retained.Contains((u, v)))
                {
                    long sum = 0; foreach (var table in tables) sum += table.Score(u, v);
                    var c = new Candidate { U = u, V = v, Score = (double)sum / Scale };
                    var priority = (c.Score, -v, -u);
                    if (queue.Count < options.Keep) Push(c);
                    else if (queue.TryPeek(out _, out var worst) && priority.CompareTo(worst) > 0)
                    { var removed = queue.Dequeue(); retained.Remove((removed.U, removed.V)); Push(c); }
                }
                if (state.PrimitiveDraws % 100000 == 0)
                { Checkpoint(); Console.WriteLine($"sample {state.PrimitiveDraws}/{options.Samples}: {timer.Elapsed.TotalSeconds:F2}s"); }
            }
        }
        finally { Checkpoint(); }
        var candidates = state.Retained.Select(c => new Candidate { U = c.U, V = c.V }).ToList();
        Hunt.ScoreStage(candidates, 0, 4093, options.Workers, token);
        candidates = Ordered(candidates, options.RefineKeep).ToList();
        Hunt.Save(Path.Combine(directory, "stage-4093.json"), candidates);
        Hunt.ScoreStage(candidates, 4093, 16381, options.Workers, token);
        candidates = Ordered(candidates, options.FinalKeep).ToList();
        foreach (var c in candidates) c.ScreeningScore = c.Score;
        // Freeze this list before the final prime interval is scored.
        Hunt.Save(Path.Combine(directory, "stage-16381.json"), candidates);
        Hunt.ScoreStage(candidates, 16381, options.PrimeBound, options.Workers, token);
        foreach (var c in candidates) c.ValidationScore = c.Score - c.ScreeningScore;
        candidates = Ordered(candidates, candidates.Count).ToList();
        var seenJ = new HashSet<string> { Family302.Curve(164518, 924945).JInvariant.ToString() };
        var output = new List<object>(); var excluded = new List<object>();
        foreach (var c in candidates)
        {
            token.ThrowIfCancellationRequested();
            var (curve, points) = Structured302.Create(c.U, c.V);
            string j = curve.JInvariant.ToString();
            if (!seenJ.Add(j)) { excluded.Add(new { c.U, c.V, reason = "duplicate_or_record_j_invariant" }); continue; }
            var cert = curve.GetRankLowerBound(points, 1009, token);
            string file = $"curve_{c.U}_{c.V}.json";
            Hunt.Save(Path.Combine(directory, file), Hunt.CurveData(curve, points,
                $"icarm302-17 generic sections at {c.U}/{c.V}; generated by bounded parameter sampler", cert.LowerBound));
            output.Add(new { id = $"{c.U}_{c.V}", u = c.U, v = c.V, file,
                score = c.Score, screening_score = c.ScreeningScore, tail_score = c.ValidationScore,
                seed_lower_bound = cert.LowerBound, j_invariant = j, published_novelty_verified = false });
        }
        Hunt.Save(Path.Combine(directory, "candidates.json"), output);
        Hunt.Save(Path.Combine(directory, "complete.json"), new { options, state.Draws, state.PrimitiveDraws,
            retained_unique_parameters = state.Retained.Length, exported_curves = output.Count, excluded,
            seconds_this_session = timer.Elapsed.TotalSeconds, score_is_not_a_rank_bound = true,
            sampling_with_replacement = true, record_j_excluded = true, other_published_curves_not_checked = true });
        Console.WriteLine($"Prepared {output.Count} different candidate j-invariants in {timer.Elapsed.TotalSeconds:F2}s.");
    }
}
