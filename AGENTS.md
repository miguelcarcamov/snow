# AGENTS.md

Instructions for AI coding agents (Cursor, Claude Code, Codex, …) working on snow.

1. **Read the conventions first:** [`.cursor/rules/`](.cursor/rules/README.md), mirrored in [`CONVENTIONS.md`](CONVENTIONS.md). Start with `general.mdc` (what snow is for) and `selfcal-domain.mdc` (self-cal invariants that must never be broken).
2. **Work from a GitHub issue.** The roadmap is milestone "v1.0 redesign" (epic #22). Stay within the issue's scope, and reference it in commits and PRs.
3. **Branch and PR:** `feature/<topic>` or `fix/<topic>` from `development`, with a PR into `development`. Never push to `development` or `master` directly.
4. **Before opening a PR:** `pre-commit run --all-files`, `pytest tests/unit/`, and `python tools/sync_conventions.py --check` if you touched `.cursor/rules/`.
5. **Science first:** if a change affects how gains are solved or applied, how the model is written, or how a round is accepted, explain the theory in the PR and add or extend the simulated end-to-end test (#32).

Claude Code-specific notes are in [`CLAUDE.md`](CLAUDE.md).
