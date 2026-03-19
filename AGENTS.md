# Repository Guidelines

## Project Structure & Module Organization
Stellerator is a single-crate Rust project. `Cargo.toml` defines the crate and dependencies, and `src/main.rs` is the current entry point. Add new modules under `src/` and declare them from `main.rs` or a future `lib.rs`. Put cross-module integration tests in `tests/` when they are needed. Generated build artifacts belong in `target/` and must stay out of version control.

## Build, Test, and Development Commands
Run commands from the repository root.

- `cargo run` builds and starts the application locally.
- `cargo build` compiles the crate without running it.
- `cargo test` runs unit and integration tests.
- `cargo fmt` formats the code with `rustfmt`.
- `cargo fmt -- --check` verifies formatting for CI or pre-PR checks.
- `cargo clippy --all-targets --all-features -- -D warnings` runs linting and treats warnings as failures.

## Coding Style & Naming Conventions
Follow standard Rust style with 4-space indentation. Use `snake_case` for functions, variables, and module files; `PascalCase` for structs, enums, and traits; and `SCREAMING_SNAKE_CASE` for constants. Keep modules small and focused. Prefer explicit public APIs and add short `///` doc comments where behavior is not obvious.

## Testing Guidelines
Use Rust’s built-in test framework. Keep unit tests in `#[cfg(test)]` modules beside the code they cover, and use `tests/` for behavior that spans modules. Name tests after the behavior they verify, such as `parses_empty_input` or `rejects_invalid_state`. Run `cargo test` before committing. New features should include tests unless there is a clear reason they are not practical yet.

## Commit & Pull Request Guidelines
The visible history is minimal, so use short, imperative commit subjects such as `Add orbit calculation helper`. Keep each commit focused on one change. Pull requests should include a brief summary, linked `bd` issue, and the commands used for validation, typically `cargo test`, `cargo fmt -- --check`, and `cargo clippy --all-targets --all-features -- -D warnings`. Include sample output when behavior changes are user-visible.

## Issue Tracking & Configuration
Track all work in `bd`; do not add Markdown TODO lists. Check available work with `bd ready --json`, claim with `bd update <id> --claim --json`, and close completed work with `bd close <id> --reason "Completed" --json`. Keep secrets out of the repository and prefer environment variables for local configuration.

<!-- BEGIN BEADS INTEGRATION v:1 profile:full hash:d4f96305 -->
## Issue Tracking with bd (beads)

**IMPORTANT**: This project uses **bd (beads)** for ALL issue tracking. Do NOT use markdown TODOs, task lists, or other tracking methods.

### Why bd?

- Dependency-aware: Track blockers and relationships between issues
- Git-friendly: Dolt-powered version control with native sync
- Agent-optimized: JSON output, ready work detection, discovered-from links
- Prevents duplicate tracking systems and confusion

### Quick Start

**Check for ready work:**

```bash
bd ready --json
```

**Create new issues:**

```bash
bd create "Issue title" --description="Detailed context" -t bug|feature|task -p 0-4 --json
bd create "Issue title" --description="What this issue is about" -p 1 --deps discovered-from:bd-123 --json
```

**Claim and update:**

```bash
bd update <id> --claim --json
bd update bd-42 --priority 1 --json
```

**Complete work:**

```bash
bd close bd-42 --reason "Completed" --json
```

### Issue Types

- `bug` - Something broken
- `feature` - New functionality
- `task` - Work item (tests, docs, refactoring)
- `epic` - Large feature with subtasks
- `chore` - Maintenance (dependencies, tooling)

### Priorities

- `0` - Critical (security, data loss, broken builds)
- `1` - High (major features, important bugs)
- `2` - Medium (default, nice-to-have)
- `3` - Low (polish, optimization)
- `4` - Backlog (future ideas)

### Workflow for AI Agents

1. **Check ready work**: `bd ready` shows unblocked issues
2. **Claim your task atomically**: `bd update <id> --claim`
3. **Work on it**: Implement, test, document
4. **Discover new work?** Create linked issue:
   - `bd create "Found bug" --description="Details about what was found" -p 1 --deps discovered-from:<parent-id>`
5. **Complete**: `bd close <id> --reason "Done"`

### Auto-Sync

bd automatically syncs via Dolt:

- Each write auto-commits to Dolt history
- Use `bd dolt push`/`bd dolt pull` for remote sync
- No manual export/import needed!

### Important Rules

- ✅ Use bd for ALL task tracking
- ✅ Always use `--json` flag for programmatic use
- ✅ Link discovered work with `discovered-from` dependencies
- ✅ Check `bd ready` before asking "what should I work on?"
- ❌ Do NOT create markdown TODO lists
- ❌ Do NOT use external issue trackers
- ❌ Do NOT duplicate tracking systems

For more details, see README.md and docs/QUICKSTART.md.

## Landing the Plane (Session Completion)

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds

<!-- END BEADS INTEGRATION -->
