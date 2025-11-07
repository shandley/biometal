# Release Instructions for v1.1.0

## ✅ Pre-Release Checklist (COMPLETED)

All automated preparation is done:

- ✅ Version bumped to 1.1.0 (Cargo.toml, pyproject.toml)
- ✅ CHANGELOG.md updated with comprehensive v1.1.0 section
- ✅ CLAUDE.md project structure updated
- ✅ BgzipWriter visibility warning fixed
- ✅ All tests passing (260 tests)
- ✅ Python wheels building successfully
- ✅ Release notes prepared (RELEASE_NOTES_v1.1.0.md)
- ✅ All commits pushed to main

---

## 📋 Manual Release Steps

### Step 1: Publish to crates.io (Rust) 🦀

```bash
# Ensure you're logged in to crates.io
cargo login

# Publish to crates.io
cargo publish

# Verify publication (takes ~1 minute to appear)
open https://crates.io/crates/biometal
```

**Expected output**:
```
Updating crates.io index
   Packaging biometal v1.1.0 (/Users/scotthandley/Code/biometal)
   Verifying biometal v1.1.0 (/Users/scotthandley/Code/biometal)
   Uploading biometal v1.1.0 (/Users/scotthandley/Code/biometal)
```

---

### Step 2: Create GitHub Release 🏷️

**This will automatically trigger PyPI publishing via GitHub Actions!**

1. **Go to GitHub Releases**:
   - URL: https://github.com/shandley/biometal/releases/new

2. **Fill in Release Details**:
   - **Tag version**: `v1.1.0`
   - **Target**: `main`
   - **Release title**: `v1.1.0 - K-mer Operations & Complexity`

3. **Copy Release Notes**:
   - Open `RELEASE_NOTES_v1.1.0.md` and copy the content
   - Paste into the description field

4. **Publish Release**:
   - Click **"Publish release"** button

---

### Step 3: Verify Automated PyPI Publishing 🐍

The GitHub Action will automatically:
- Build wheels for macOS ARM, macOS x86_64, Linux x86_64
- Build source distribution
- Publish to PyPI using trusted publishing (no token needed!)

**Monitor the GitHub Action**:
1. Go to: https://github.com/shandley/biometal/actions
2. Look for workflow: "Publish to PyPI"
3. Watch progress (takes ~10-15 minutes)

**Verify PyPI publication** (after GitHub Action completes):
```bash
# Check PyPI (takes ~5 minutes after workflow completes)
open https://pypi.org/project/biometal-rs/

# Test installation
pip install biometal-rs==1.1.0

# Verify version
python -c "import biometal; print(biometal.__version__)"
# Expected: 1.1.0
```

---

## 🔍 Verification Checklist

After publishing, verify:

### Crates.io ✅
- [ ] https://crates.io/crates/biometal shows v1.1.0
- [ ] Documentation at https://docs.rs/biometal/1.1.0/ is live
- [ ] `cargo search biometal` shows v1.1.0

### PyPI ✅
- [ ] https://pypi.org/project/biometal-rs/ shows v1.1.0
- [ ] Wheels available for:
  - [ ] macOS ARM (aarch64)
  - [ ] macOS x86_64
  - [ ] Linux x86_64
- [ ] Source distribution (.tar.gz) available
- [ ] `pip search biometal-rs` shows v1.1.0 (if pip search is enabled)

### GitHub ✅
- [ ] Release v1.1.0 visible at https://github.com/shandley/biometal/releases
- [ ] Release notes properly formatted
- [ ] GitHub Action "Publish to PyPI" completed successfully

---

## 🎉 Post-Release

### Announce Release

**Options**:
1. **Twitter/X**: Share release with #bioinformatics #rust #python
2. **Reddit**: Post to r/bioinformatics, r/rust
3. **Discord**: Share in Rust and bioinformatics communities
4. **Email**: Notify any early adopters or collaborators

**Sample announcement**:
```
🧬 biometal v1.1.0 released!

New features:
- K-mer operations with evidence-based design (Entry 034)
- Complexity scoring (Shannon entropy)
- Python bindings for all k-mer functions
- Grade A+ code quality

Evidence-based: K-mers are data-structure-bound, not compute-bound.
NEON provides no benefit - scalar is optimal!

Rust: https://crates.io/crates/biometal
Python: https://pypi.org/project/biometal-rs/
```

### Monitor

- Watch for GitHub issues
- Monitor PyPI download stats
- Check crates.io download stats

---

## 🐛 Troubleshooting

### Cargo publish fails

**Issue**: "error: failed to verify..."
```bash
# Check for uncommitted changes
git status

# Re-run tests
cargo test --lib

# Try dry-run first
cargo publish --dry-run
```

### GitHub Action fails

**Common issues**:
1. **PyPI trusted publishing not configured**:
   - Go to https://pypi.org/manage/account/publishing/
   - Add trusted publisher for shandley/biometal

2. **Build failures**:
   - Check workflow logs
   - Verify Cargo.toml and pyproject.toml versions match

3. **Permission errors**:
   - Verify PyPI environment is configured in GitHub repo settings

### PyPI wheels missing

**Issue**: Some wheels didn't upload

```bash
# Manually build and upload if needed
maturin build --release
maturin publish

# Or trigger GitHub Action manually
# Go to Actions → Publish to PyPI → Run workflow
```

---

## 📊 Expected Timeline

- **Step 1 (crates.io)**: ~2-3 minutes
- **Step 2 (GitHub release)**: ~1 minute
- **Step 3 (GitHub Action → PyPI)**: ~10-15 minutes
  - Build wheels: ~8-10 minutes
  - Upload to PyPI: ~1-2 minutes
  - PyPI indexing: ~3-5 minutes

**Total**: ~15-20 minutes from start to finish

---

## ✅ Success Criteria

Release is successful when:
1. ✅ `cargo search biometal` shows v1.1.0
2. ✅ https://crates.io/crates/biometal shows v1.1.0
3. ✅ https://pypi.org/project/biometal-rs/ shows v1.1.0
4. ✅ `pip install biometal-rs==1.1.0` works
5. ✅ GitHub release v1.1.0 is published
6. ✅ All GitHub Actions workflows pass

---

**Good luck with the release! 🚀**

If you encounter any issues, check the troubleshooting section or ping me for help.
