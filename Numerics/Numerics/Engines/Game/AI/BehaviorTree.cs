using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Game.AI;

/// <summary>
/// Lightweight behavior tree executor for hybrid AI.
/// Supports Sequence, Selector, Condition, Action, and Decorator nodes.
///
/// Designed for combining ML-driven decisions with scripted fallback logic.
/// For example:
///   Selector(
///     Sequence(IsInCombatRange, MLDogfightAction),
///     Sequence(IsLowFuel, ReturnToBase),
///     PatrolWaypoints)
///
/// Execution is tick-based: call Tick(context) each frame.
/// </summary>
public class BehaviorTree
{
    private readonly BehaviorNode _root;

    /// <summary>Name of this behavior tree.</summary>
    public string Name { get; }

    /// <summary>
    /// Creates a behavior tree with the given root node.
    /// </summary>
    public BehaviorTree(string name, BehaviorNode root)
    {
        Name = name;
        _root = root ?? throw new ArgumentNullException(nameof(root));
    }

    /// <summary>
    /// Tick the behavior tree once with the given context.
    /// Returns the result of the root node's execution.
    /// </summary>
    public NodeStatus Tick(BehaviorContext context)
    {
        return _root.Execute(context);
    }

    /// <summary>
    /// Reset all node state (e.g. running nodes).
    /// </summary>
    public void Reset()
    {
        _root.Reset();
    }
}

/// <summary>
/// Result of a behavior tree node execution.
/// </summary>
public enum NodeStatus
{
    /// <summary>Node completed successfully.</summary>
    Success,

    /// <summary>Node failed.</summary>
    Failure,

    /// <summary>Node is still running (will be resumed next tick).</summary>
    Running
}

/// <summary>
/// Shared context passed to all behavior tree nodes during execution.
/// Contains arbitrary key-value pairs for game state.
/// </summary>
public class BehaviorContext
{
    private readonly Dictionary<string, object> _data = new();

    /// <summary>Set a value in the context.</summary>
    public void Set<T>(string key, T value) => _data[key] = value;

    /// <summary>Get a value from the context.</summary>
    public T Get<T>(string key) => _data.TryGetValue(key, out var val) ? (T)val : default;

    /// <summary>Check if a key exists.</summary>
    public bool Has(string key) => _data.ContainsKey(key);

    /// <summary>Get the raw data dictionary.</summary>
    public Dictionary<string, object> Data => _data;
}

/// <summary>
/// Base class for all behavior tree nodes.
/// </summary>
public abstract class BehaviorNode
{
    /// <summary>Optional node name for debugging.</summary>
    public string Name { get; set; }

    /// <summary>Execute this node and return its status.</summary>
    public abstract NodeStatus Execute(BehaviorContext context);

    /// <summary>Reset any running state.</summary>
    public virtual void Reset() { }
}

/// <summary>
/// Sequence node: executes children left-to-right.
/// Returns Success if all children succeed.
/// Returns Failure as soon as any child fails.
/// Returns Running if a child is still running (resumes there next tick).
/// </summary>
public class SequenceNode : BehaviorNode
{
    private readonly List<BehaviorNode> _children;
    private int _currentIndex;

    public SequenceNode(string name, params BehaviorNode[] children)
    {
        Name = name;
        _children = new List<BehaviorNode>(children);
    }

    public override NodeStatus Execute(BehaviorContext context)
    {
        while (_currentIndex < _children.Count)
        {
            var status = _children[_currentIndex].Execute(context);
            if (status == NodeStatus.Failure)
            {
                _currentIndex = 0;
                return NodeStatus.Failure;
            }
            if (status == NodeStatus.Running)
                return NodeStatus.Running;

            _currentIndex++;
        }

        _currentIndex = 0;
        return NodeStatus.Success;
    }

    public override void Reset()
    {
        _currentIndex = 0;
        foreach (var child in _children) child.Reset();
    }
}

/// <summary>
/// Selector node: executes children left-to-right.
/// Returns Success as soon as any child succeeds.
/// Returns Failure if all children fail.
/// Returns Running if a child is still running.
/// </summary>
public class SelectorNode : BehaviorNode
{
    private readonly List<BehaviorNode> _children;
    private int _currentIndex;

    public SelectorNode(string name, params BehaviorNode[] children)
    {
        Name = name;
        _children = new List<BehaviorNode>(children);
    }

    public override NodeStatus Execute(BehaviorContext context)
    {
        while (_currentIndex < _children.Count)
        {
            var status = _children[_currentIndex].Execute(context);
            if (status == NodeStatus.Success)
            {
                _currentIndex = 0;
                return NodeStatus.Success;
            }
            if (status == NodeStatus.Running)
                return NodeStatus.Running;

            _currentIndex++;
        }

        _currentIndex = 0;
        return NodeStatus.Failure;
    }

    public override void Reset()
    {
        _currentIndex = 0;
        foreach (var child in _children) child.Reset();
    }
}

/// <summary>
/// Condition node: evaluates a predicate and returns Success or Failure.
/// Never returns Running.
/// </summary>
public class ConditionNode : BehaviorNode
{
    private readonly Func<BehaviorContext, bool> _predicate;

    public ConditionNode(string name, Func<BehaviorContext, bool> predicate)
    {
        Name = name;
        _predicate = predicate;
    }

    public override NodeStatus Execute(BehaviorContext context)
    {
        return _predicate(context) ? NodeStatus.Success : NodeStatus.Failure;
    }
}

/// <summary>
/// Action node: executes a function and returns its status.
/// The function should return Success, Failure, or Running.
/// </summary>
public class ActionNode : BehaviorNode
{
    private readonly Func<BehaviorContext, NodeStatus> _action;

    public ActionNode(string name, Func<BehaviorContext, NodeStatus> action)
    {
        Name = name;
        _action = action;
    }

    public override NodeStatus Execute(BehaviorContext context)
    {
        return _action(context);
    }
}

/// <summary>
/// Inverter decorator: inverts the result of its child.
/// Success → Failure, Failure → Success, Running → Running.
/// </summary>
public class InverterNode : BehaviorNode
{
    private readonly BehaviorNode _child;

    public InverterNode(string name, BehaviorNode child)
    {
        Name = name;
        _child = child;
    }

    public override NodeStatus Execute(BehaviorContext context)
    {
        var status = _child.Execute(context);
        return status switch
        {
            NodeStatus.Success => NodeStatus.Failure,
            NodeStatus.Failure => NodeStatus.Success,
            _ => NodeStatus.Running
        };
    }

    public override void Reset() => _child.Reset();
}

/// <summary>
/// Repeater decorator: repeats its child N times (or until failure).
/// Returns Success after N successful executions.
/// </summary>
public class RepeaterNode : BehaviorNode
{
    private readonly BehaviorNode _child;
    private readonly int _count;
    private int _current;

    /// <param name="name">Node name.</param>
    /// <param name="child">Child node to repeat.</param>
    /// <param name="count">Number of repetitions. 0 = infinite.</param>
    public RepeaterNode(string name, BehaviorNode child, int count = 0)
    {
        Name = name;
        _child = child;
        _count = count;
    }

    public override NodeStatus Execute(BehaviorContext context)
    {
        var status = _child.Execute(context);
        if (status == NodeStatus.Running) return NodeStatus.Running;
        if (status == NodeStatus.Failure) { _current = 0; return NodeStatus.Failure; }

        _current++;
        if (_count > 0 && _current >= _count)
        {
            _current = 0;
            return NodeStatus.Success;
        }

        return NodeStatus.Running;
    }

    public override void Reset()
    {
        _current = 0;
        _child.Reset();
    }
}
