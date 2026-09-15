"""A bounded mixture of short and more widely supported known-point vectors."""
import random

from bounded_anchor_pool import bounded_vectors


def diverse_vectors(gram, limit, count):
    n=len(gram)
    if n<3:return bounded_vectors(gram,limit,count)
    rows,short_stats=bounded_vectors(gram,max(n,limit//2),(count+1)//2)
    chosen={v for _,v in rows};sample={};rng=random.Random(20260914)
    remaining=max(0,limit-short_stats['vectors_considered'])
    # These support sizes are fixed in advance, independent of any target point.
    supports=list(range(2,min(8,n)+1));attempts=0
    while len(sample)<remaining and attempts<4*max(1,remaining):
        support=supports[attempts%len(supports)];attempts+=1
        v=[0]*n
        for i in rng.sample(range(n),support):v[i]=rng.choice([-1,1])
        if next(x for x in v if x)<0:v=[-x for x in v]
        v=tuple(v)
        if v in sample or v in chosen:continue
        indices=[i for i in range(n) if v[i]]
        score=sum(v[i]*gram[i][j]*v[j] for i in indices for j in indices)
        sample[v]=(score,support)
    groups={k:[] for k in supports}
    for v,(score,k) in sample.items():groups[k].append((score,v))
    for group in groups.values():group.sort(reverse=True)
    while len(rows)<count and any(groups.values()):
        for k in supports:
            if groups[k] and len(rows)<count:rows.append(groups[k].pop())
    rows.sort()
    return rows,{'vectors_considered':short_stats['vectors_considered']+len(sample),
        'vector_limit':limit,'short_vectors_considered':short_stats['vectors_considered'],
        'random_vectors_considered':len(sample),'random_attempts':attempts,
        'random_attempt_limit':4*max(1,remaining),'random_supports':supports,
        'coefficient_bound':3,'l1_bound':12,'selection':'half short, half support strata',
        'selected_support_counts':{str(k):sum(sum(x!=0 for x in v)==k for _,v in rows)
                                   for k in range(1,n+1)}}
