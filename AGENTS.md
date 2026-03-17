# Repository Guidelines

## Project Structure & Module Organization
`Stellerator` is currently a single-crate Rust project. The crate manifest lives in `Cargo.toml`, and the application entry point is `src/main.rs`. Add new modules under `src/` and declare them from `main.rs` or `lib.rs` if the project is later split into a library and binary. Rust build artifacts are generated in `target/`; do not commit that directory.

## Build, Test, and Development Commands
- `cargo run`: build and launch the application locally.
- `cargo build`: compile the crate without running it.
- `cargo test`: run unit and integration tests.
- `cargo fmt`: format the codebase with `rustfmt`.
- `cargo fmt -- --check`: verify formatting in CI or before opening a PR.
- `cargo clippy --all-targets --all-features -- -D warnings`: run lint checks and fail on warnings.

Use these commands from the repository root.

## Coding Style & Naming Conventions
Follow standard Rust style: 4-space indentation, one item per line when lists grow, and `snake_case` for functions, variables, and module files. Use `PascalCase` for types and traits and `SCREAMING_SNAKE_CASE` for constants. Prefer small modules with focused responsibilities. Keep public APIs explicit and document non-obvious behavior with concise `///` doc comments.

## Testing Guidelines
Use Rust’s built-in test framework. Place unit tests in `#[cfg(test)]` modules beside the code they cover, and add integration tests under `tests/` when behavior crosses module boundaries. Name tests for the behavior they verify, for example `parses_empty_input` or `rejects_invalid_state`. Run `cargo test` before every commit. New features should ship with tests or a clear justification for why tests are not practical yet.

## Commit & Pull Request Guidelines
This repository does not have commit history yet, so establish a clean convention now: use short, imperative commit subjects such as `Add orbit calculation helper`. Keep commits focused and logically grouped. Pull requests should include a brief summary, test notes (`cargo test`, `cargo clippy`, `cargo fmt -- --check`), and any screenshots or sample output if behavior changes are user-visible.

## Configuration & Hygiene
Keep secrets out of the repository. Prefer environment variables for local configuration, and document any new setup requirements in `README.md` or this guide when the project grows.

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
