"""
Entry point for the GeneRNBI agentic assistant (Biomni A1 orchestrator).
"""

import sys
from pathlib import Path

# Allow running from src/ directly
sys.path.insert(0, str(Path(__file__).parent))

from agent import initialize_agent


def interactive_mode():
    """Run the agent in interactive mode."""
    print("\n" + "=" * 70)
    print("GeneRNBI Assistant — Powered by Biomni A1")
    print("=" * 70)
    print("\nI can help you with:")
    print("  • Understanding the benchmark pipeline and metrics")
    print("  • Running evaluations on datasets (Norman, Replogle, etc.)")
    print("  • Integrating new GRN inference methods")
    print("  • Debugging errors and configuration issues")
    print("\nType 'quit' or 'exit' to end the session.")
    print("=" * 70 + "\n")

    agent = initialize_agent()

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
            _log, _answer = agent.go(user_input)
            print(f"\n{'=' * 70}")
        except Exception as e:
            print(f"\n❌ Error: {e}\n")


def single_query(question: str):
    """Run a single query."""
    print(f"\nQuestion: {question}\n{'=' * 70}")
    agent = initialize_agent()
    _log, _answer = agent.go(question)
    print(f"\n{'=' * 70}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        single_query(" ".join(sys.argv[1:]))
    else:
        interactive_mode()
