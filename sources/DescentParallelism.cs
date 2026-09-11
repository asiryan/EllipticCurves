using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.ExceptionServices;
using System.Threading;
using System.Threading.Tasks;

namespace EllipticCurves
{
    internal static class DescentParallelism
    {
        internal static void ForEach<T>(DescentBudget budget,
            Func<DescentBudget, IEnumerable<T>> items, Action<T, DescentBudget> action)
        {
            if (budget.Options.MaxDegreeOfParallelism == 1)
            {
                foreach (var item in items(budget)) action(item, budget);
                return;
            }

            using (var stop = CancellationTokenSource.CreateLinkedTokenSource(budget.Token))
            {
                var workerBudget = budget.WithToken(stop.Token);
                ExceptionDispatchInfo failure = null;
                ExceptionDispatchInfo unexpectedFailure = null;
                // Stop peers even when the failure comes from MoveNext, rather than a body.
                void RecordFailure(Exception error)
                {
                    // A simultaneous arithmetic/programming error must never be hidden
                    // by a work-limit exception and returned as an ordinary partial bound.
                    if (!(error is DescentLimitException) && !(error is OperationCanceledException && stop.IsCancellationRequested))
                        Interlocked.CompareExchange(ref unexpectedFailure, ExceptionDispatchInfo.Capture(error), null);
                    Interlocked.CompareExchange(ref failure, ExceptionDispatchInfo.Capture(error), null);
                    stop.Cancel();
                }
                IEnumerable<T> ReadItems()
                {
                    using (var iterator = items(workerBudget).GetEnumerator())
                    {
                        while (true)
                        {
                            bool hasNext;
                            try { workerBudget.Token.ThrowIfCancellationRequested(); hasNext = iterator.MoveNext(); }
                            catch (Exception error) { RecordFailure(error); throw; }
                            if (!hasNext) yield break;
                            yield return iterator.Current;
                        }
                    }
                }

                try
                {
                    // NoBuffering bounds speculative enumeration and distributes uneven rows
                    // dynamically, without materializing the potentially enormous search region.
                    Parallel.ForEach(Partitioner.Create(ReadItems(), EnumerablePartitionerOptions.NoBuffering),
                        new ParallelOptions { MaxDegreeOfParallelism = budget.Options.MaxDegreeOfParallelism,
                            CancellationToken = workerBudget.Token }, item =>
                        {
                            try { workerBudget.Token.ThrowIfCancellationRequested(); action(item, workerBudget); }
                            catch (Exception error) { RecordFailure(error); throw; }
                        });
                }
                catch (Exception error) when (error is AggregateException || error is OperationCanceledException)
                {
                    // Parallel.ForEach joins every worker before throwing. Only then may the
                    // caller read the partial proof. Preserve library cancellation/limit types.
                    if (unexpectedFailure != null) unexpectedFailure.Throw();
                    budget.Token.ThrowIfCancellationRequested();
                    if (failure != null) failure.Throw();
                    if (error is AggregateException aggregate && aggregate.Flatten().InnerExceptions.Count == 1)
                        ExceptionDispatchInfo.Capture(aggregate.Flatten().InnerExceptions.Single()).Throw();
                    throw;
                }
                budget.Token.ThrowIfCancellationRequested();
            }
        }
    }
}
