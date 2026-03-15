"""
Entry point for the GeneRNBI agentic assistant.

Usage:
    python run.py [--agent chat|a1] ["question"]

Agent modes:
    chat  (default) — lightweight chat completion; best for Q&A, debugging,
                       and knowledge lookup. Uses know-how docs + search_viash tool.
    a1              — Biomni A1 task-execution agent; use when you need to run
                       commands, execute Viash pipelines, or modify files.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))


def _make_agent(mode: str):
    if mode == "a1":
        from agent import initialize_agent
        return initialize_agent(), "a1"
    else:
        from chat_agent import initialize_chat_agent
        return initialize_chat_agent(), "chat"


def _run_chat(agent, question: str):
    answer = agent.ask(question)
    print(f"\n{answer}\n")


def _run_a1(agent, question: str):
    _log, _answer = agent.go(question)


def interactive_mode(mode: str):
    label = "Chat" if mode == "chat" else "Biomni A1"
    print("\n" + "=" * 70)
    print(f"GeneRNBI Assistant — {label} mode")
    print("=" * 70)
    print("\nI can help you with:")
    print("  • Understanding the benchmark pipeline and metrics")
    print("  • Debugging errors and configuration issues")
    print("  • Integrating new GRN inference methods")
    print("  • Running evaluations on datasets")
    if mode == "chat":
        print("\nTip: use --agent a1 to switch to task-execution mode.")
    print("\nType 'quit' or 'exit' to end the session.")
    print("=" * 70 + "\n")

    agent, resolved_mode = _make_agent(mode)

    while True:
        try:
            user_input = input("\nYou: ").strip()
        except (KeyboardInterrupt, EOFError):
            print("\n\nExiting...")
            break

        if not user_input:
            continue
        if user_input.lower() in ["quit", "exit"]:
            print("\nGoodbye!")
            break

        try:
            print(f"\n{'=' * 70}")
            if resolved_mode == "chat":
                _run_chat(agent, user_input)
            else:
                _run_a1(agent, user_input)
            print(f"{'=' * 70}")
        except Exception as e:
            print(f"\n❌ Error: {e}\n")


def single_query(question: str, mode: str):
    print(f"\nQuestion: {question}\n{'=' * 70}")
    agent, resolved_mode = _make_agent(mode)
    if resolved_mode == "chat":
        _run_chat(agent, question)
    else:
        _run_a1(agent, question)
    print(f"{'=' * 70}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--agent", choices=["chat", "a1"], default="chat")
    args, remaining = parser.parse_known_args()

    question = " ".join(remaining).strip() if remaining else None

    if question:
        single_query(question, args.agent)
    else:
        interactive_mode(args.agent)
