using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Common
{
    /// <summary>
    /// Lightweight publish/subscribe event bus for engine events.
    /// Events are keyed by type; subscribers receive instances of the event type.
    /// </summary>
    public class EventBus
    {
        private readonly Dictionary<Type, List<Delegate>> _handlers = new();

        /// <summary>Subscribe to events of type <typeparamref name="T"/>.</summary>
        public void Subscribe<T>(Action<T> handler)
        {
            var key = typeof(T);
            if (!_handlers.TryGetValue(key, out var list))
            {
                list = new List<Delegate>();
                _handlers[key] = list;
            }
            list.Add(handler);
        }

        /// <summary>Unsubscribe a previously registered handler.</summary>
        public void Unsubscribe<T>(Action<T> handler)
        {
            var key = typeof(T);
            if (_handlers.TryGetValue(key, out var list))
                list.Remove(handler);
        }

        /// <summary>Publish an event to all subscribers of type <typeparamref name="T"/>.</summary>
        public void Publish<T>(T evt)
        {
            var key = typeof(T);
            if (_handlers.TryGetValue(key, out var list))
            {
                foreach (var handler in list)
                    ((Action<T>)handler)(evt);
            }
        }

        /// <summary>Remove all subscribers.</summary>
        public void Clear() => _handlers.Clear();
    }
}
