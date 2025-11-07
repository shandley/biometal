# Documentation Update Plan - Post v1.2.0 Release

**Date**: November 6, 2025
**Purpose**: Update all documentation to reflect v1.2.0 release status
**Scope**: Version references, status updates, roadmap, feature completeness

---

## 1. CLAUDE.md Updates

### Current Issues
- Header shows "v1.0.0 - Production Release"
- Project Status section mentions v1.0.0 as current
- Test counts outdated (279 → 347 with v1.2.0)
- Current Work section shows "v1.1.0 - K-mer Operations"
- Missing v1.2.0 Python bindings in status
- "Next focus" section outdated

### Required Changes
- [ ] Update header: `v1.0.0` → `v1.2.0`
- [ ] Update status section with v1.2.0 release info
- [ ] Update test counts: 279 → 347 (260 library + 87 doc)
- [ ] Add Python module structure (src/python/)
- [ ] Update "Current Work" → "Recent Releases"
- [ ] Update "Next focus" section
- [ ] Update "Last Updated" timestamp

---

## 2. README.md Updates

### Current Issues
- Roadmap shows "v1.0.0 (Released)" as latest
- Missing v1.1.0 and v1.2.0 in roadmap
- Status footer shows "v1.0.0"
- Test counts outdated (121 → 347)
- Quick Start version examples show "1.0"

### Required Changes
- [ ] Update Quick Start: `biometal = "1.0"` → `biometal = "1.2"`
- [ ] Update Roadmap section:
  - v1.0.0 ✅ Released
  - v1.1.0 ✅ Released (K-mer Operations)
  - v1.2.0 ✅ Released (Python Phase 4 Bindings)
  - Future Considerations (unchanged)
- [ ] Update status footer: v1.0.0 → v1.2.0
- [ ] Update test counts: 121 → 347
- [ ] Update grade: "A+" → "A" (current)

---

## 3. Project Structure Documentation

### Required Changes
- [ ] Add src/python/ module structure in CLAUDE.md:
  ```
  ├── src/python/          # Python bindings (PyO3 0.27)
  │   ├── mod.rs           # Module registration
  │   ├── records.rs       # PyFastqRecord, PyFastaRecord
  │   ├── streams.rs       # PyFastqStream, PyFastaStream
  │   ├── operations.rs    # Core operations (GC, base counting)
  │   ├── kmers.rs         # K-mer operations (v1.1.0)
  │   ├── sequence.rs      # Sequence ops (v1.2.0)
  │   ├── record_ops.rs    # Record manipulation (v1.2.0)
  │   ├── trimming.rs      # Trimming operations (v1.2.0)
  │   └── masking.rs       # Masking operations (v1.2.0)
  ```

---

## 4. Version Reference Updates

### Files to Check
- [ ] Cargo.toml (already updated to 1.2.0)
- [ ] pyproject.toml (already updated to 1.2.0)
- [ ] lib.rs module docstring
- [ ] src/python/mod.rs docstring

---

## 5. Archive Old Release Files

### Files to Move to `archive/releases/`
- [ ] RELEASE_INSTRUCTIONS_v1.1.0.md
- [ ] RELEASE_NOTES_v1.0.0.md (if exists)
- [ ] RELEASE_NOTES_v1.1.0.md (if exists)
- [ ] PYPI_SETUP_CHECKLIST.md (obsolete - using GitHub Actions)
- [ ] PYPI_RELEASE_INSTRUCTIONS.md (obsolete)
- [ ] QUICKSTART.md (duplicate of README content)

Keep current:
- RELEASE_INSTRUCTIONS_v1.2.0.md (latest)
- RELEASE_NOTES_v1.2.0.md (latest)

---

## 6. Session Restart Checklist (CLAUDE.md)

### Update Quick Context
```markdown
### Quick Context
1. **Project**: biometal v1.2.0 - ARM-native bioinformatics library
2. **Status**: v1.2.0 released (Nov 6, 2025) - Python Phase 4 bindings
3. **Tests**: 347 passing (260 library + 87 doc)
4. **Grade**: A (rust-code-quality-reviewer compatible)
5. **Philosophy**: Evidence-based optimization (1,357 experiments, 40,710 measurements)
```

### Update Current Work
```markdown
### Recent Releases
- ✅ v1.0.0 (Nov 5, 2025): Core library, FASTQ/FASTA streaming, NEON ops, network streaming
- ✅ v1.1.0 (Nov 6, 2025): K-mer operations & complexity scoring
- ✅ v1.2.0 (Nov 6, 2025): Python bindings for Phase 4 sequence operations (20 functions)

### Known State
- **Package naming**: `biometal-rs` on PyPI, `biometal` on crates.io/GitHub
- **Latest version**: v1.2.0 (both platforms)
- **Python API**: 40+ functions (core ops + k-mers + Phase 4)
- **Test coverage**: 347 tests (260 library + 87 doc)

### Next Priorities (TBD)
- Community feedback on v1.2.0
- Python examples & tutorials
- BAM/SAM format support (evidence-based evaluation needed)
- Performance benchmarking vs existing tools
```

---

## 7. Documentation Status Report

Create `DOCUMENTATION_STATUS.md` with:
- Current version: v1.2.0
- Documentation coverage audit
- Known gaps
- Future documentation needs
- Link to all release notes

---

## Implementation Order

1. **CLAUDE.md** (Development guide - highest priority)
2. **README.md** (User-facing - high visibility)
3. **Archive old releases** (Cleanup)
4. **Create DOCUMENTATION_STATUS.md** (Audit report)
5. **Commit and push**

---

## Success Criteria

- [ ] All version references show v1.2.0 as current
- [ ] Roadmap accurately reflects three releases
- [ ] Test counts match actual (347)
- [ ] Python module structure documented
- [ ] Old release files archived
- [ ] Documentation status report created
- [ ] No outdated "Current Work" references
- [ ] Clean git status after commit

---

## Time Estimate

- Analysis: 10 minutes (done)
- CLAUDE.md updates: 15 minutes
- README.md updates: 10 minutes
- Archive cleanup: 5 minutes
- Status report: 10 minutes
- Testing & commit: 10 minutes

**Total**: ~60 minutes
