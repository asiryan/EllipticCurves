using System;
using System.Threading;

namespace EllipticCurves
{
    // Retain only completed values. Unlike Lazy<T>, a cancelled or failed attempt
    // must not prevent a later caller from retrying with its own token.
    internal sealed class CachedComputation<T> where T : class
    {
        private readonly SemaphoreSlim gate = new SemaphoreSlim(1, 1);
        private T value;

        internal T Get(Func<T> compute, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            var completed = Volatile.Read(ref value);
            if (completed != null) return completed;
            gate.Wait(token);
            try
            {
                token.ThrowIfCancellationRequested();
                completed = value;
                if (completed == null)
                {
                    completed = compute();
                    token.ThrowIfCancellationRequested();
                    Volatile.Write(ref value, completed);
                }
                return completed;
            }
            finally { gate.Release(); }
        }
    }
}
