# PyPI Publishing Setup Checklist

Complete guide for publishing biometal v1.0.0 to PyPI. Follow each step in order.

**Time Estimate**: 30 minutes (15 min test + 15 min production)

---

## Phase 1: TestPyPI Testing (15 minutes)

### ☐ Step 1: Create TestPyPI Account (2 minutes)

1. Open browser: https://test.pypi.org/account/register/
2. Fill in:
   - Username: `[your choice]`
   - Email: `[your email]`
   - Password: `[secure password]`
3. Click "Sign up"
4. Check email for verification link
5. Click verification link

**Checkpoint**: You can log in to test.pypi.org ✓

---

### ☐ Step 2: Configure TestPyPI Trusted Publishing (3 minutes)

1. Log in to TestPyPI: https://test.pypi.org
2. Go to: https://test.pypi.org/manage/account/publishing/
3. Click "Add a new pending publisher"
4. Fill in the form:
   ```
   PyPI Project Name: biometal-rs
   Owner: shandley
   Repository name: biometal
   Workflow name: test-pypi-publish.yml
   Environment name: testpypi
   ```
5. Click "Add"

**What this does**: Allows GitHub Actions to publish to TestPyPI without API tokens (secure OIDC).

**Checkpoint**: You see "Pending publisher added" message ✓

---

### ☐ Step 3: Create GitHub `testpypi` Environment (2 minutes)

1. Open browser: https://github.com/shandley/biometal/settings/environments
   - (If you get 404, you need admin access to the repo)
2. Click "New environment" button
3. Name: `testpypi` (exactly, case-sensitive)
4. (Optional) Add protection rules:
   - [ ] Required reviewers
   - [ ] Wait timer
5. Click "Configure environment"

**What this does**: GitHub environment for deployment, enhances security.

**Checkpoint**: You see `testpypi` in environments list ✓

---

### ☐ Step 4: Trigger Test Workflow (1 minute)

1. Go to: https://github.com/shandley/biometal/actions
2. Click "Test Publish to TestPyPI" in left sidebar
3. Click "Run workflow" button (top right)
4. Branch: `main` (should be selected)
5. Click green "Run workflow" button

**What this does**: Starts build and publish to TestPyPI.

**Checkpoint**: You see workflow running (yellow dot) ✓

---

### ☐ Step 5: Monitor Workflow (3-5 minutes)

1. Stay on Actions page
2. Click on the running workflow (e.g., "Test Publish to TestPyPI #1")
3. Watch jobs complete:
   - ⏳ Build wheels (macos-aarch64) - ~2 min
   - ⏳ Build wheels (linux-x86_64) - ~2 min
   - ⏳ Build source distribution - ~30s
   - ⏳ Test Publish to TestPyPI - ~30s

**Expected**: All jobs complete with green ✓

**If any job fails**:
- Click on the failed job
- Read error message
- Common issues:
  - Trusted publishing not configured correctly
  - Environment name mismatch
  - First-time setup delay (wait 2 min, retry)

**Checkpoint**: All 4 jobs green ✓

---

### ☐ Step 6: Verify on TestPyPI (2 minutes)

1. Go to: https://test.pypi.org/project/biometal-rs/
2. Check:
   - [ ] Package exists
   - [ ] Version is 1.0.0
   - [ ] Download files: 3 files
     - `biometal-1.0.0-cp311-cp311-macosx_11_0_arm64.whl` (~1 MB)
     - `biometal-1.0.0-cp311-cp311-manylinux_2_17_x86_64.whl` (~0.9 MB)
     - `biometal-1.0.0.tar.gz` (~332 KB)
   - [ ] Project description shows (from README.md)

**Checkpoint**: Package visible on TestPyPI ✓

---

### ☐ Step 7: Test Installation (2 minutes)

Open terminal:

```bash
# Create test environment
python3 -m venv test_env
source test_env/bin/activate  # macOS/Linux
# OR: test_env\Scripts\activate  # Windows

# Install from TestPyPI
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple/ \
            biometal==1.0.0

# Verify version
python -c "import biometal; print(f'Version: {biometal.__version__}')"
```

**Expected output**:
```
Successfully installed biometal-1.0.0
Version: 1.0.0
```

**Checkpoint**: Installation successful ✓

---

### ☐ Step 8: Test Functionality (2 minutes)

```bash
# Test basic operations (still in test_env)
python << 'EOF'
import biometal

# Create test FASTQ data
import gzip
with gzip.open('test.fq.gz', 'wt') as f:
    f.write('@read1\n')
    f.write('ATGCATGC\n')
    f.write('+\n')
    f.write('IIIIIIII\n')

# Test streaming
stream = biometal.FastqStream.from_path('test.fq.gz')
record = next(stream)

# Test operations
gc = biometal.gc_content(bytes(record.sequence))
counts = biometal.count_bases(bytes(record.sequence))
mean_q = biometal.mean_quality(bytes(record.quality))

print(f'✓ Streaming works: {record.id}')
print(f'✓ GC content: {gc:.2%}')
print(f'✓ Base counts: {counts}')
print(f'✓ Mean quality: {mean_q:.1f}')
print('\n✅ ALL TESTS PASSED')
EOF

# Cleanup
rm test.fq.gz
deactivate
rm -rf test_env
```

**Expected output**:
```
✓ Streaming works: read1
✓ GC content: 50.00%
✓ Base counts: {'A': 2, 'C': 2, 'G': 2, 'T': 2}
✓ Mean quality: 40.0

✅ ALL TESTS PASSED
```

**Checkpoint**: All operations work ✓

---

### 🎉 Phase 1 Complete!

**If all checkpoints passed**: TestPyPI publishing works correctly!

**Next**: Proceed to Phase 2 (Production PyPI)

---

## Phase 2: Production PyPI (15 minutes)

### ☐ Step 9: Create Production PyPI Account (2 minutes)

1. Open browser: https://pypi.org/account/register/
2. Fill in (can use same credentials as TestPyPI):
   - Username: `[your choice]`
   - Email: `[your email]`
   - Password: `[secure password]`
3. Click "Sign up"
4. Check email for verification link
5. Click verification link
6. **Enable 2FA** (REQUIRED for PyPI maintainers):
   - Go to: https://pypi.org/manage/account/
   - Click "Add 2FA with authentication application"
   - Scan QR code with authenticator app
   - Enter code

**Checkpoint**: You can log in to pypi.org with 2FA ✓

---

### ☐ Step 10: Configure Production PyPI Trusted Publishing (3 minutes)

1. Log in to PyPI: https://pypi.org
2. Go to: https://pypi.org/manage/account/publishing/
3. Click "Add a new pending publisher"
4. Fill in the form:
   ```
   PyPI Project Name: biometal-rs
   Owner: shandley
   Repository name: biometal
   Workflow name: publish-pypi.yml
   Environment name: pypi
   ```
5. Click "Add"

**Note**: Workflow name is different (`publish-pypi.yml`, not `test-pypi-publish.yml`)

**Checkpoint**: You see "Pending publisher added" for biometal ✓

---

### ☐ Step 11: Create GitHub `pypi` Environment (2 minutes)

1. Go to: https://github.com/shandley/biometal/settings/environments
2. Click "New environment" button
3. Name: `pypi` (exactly, case-sensitive)
4. **Recommended**: Add protection rules:
   - ✓ Required reviewers: Add yourself
   - ✓ Wait timer: 0 minutes (or 5 if you want review time)
5. Click "Configure environment"

**Why protection rules**: Prevents accidental publishes, gives you review gate.

**Checkpoint**: You see `pypi` in environments list ✓

---

### ☐ Step 12: Verify Git Tag Exists (1 minute)

```bash
cd /path/to/biometal

# Check if v1.0.0 tag exists
git tag -l v1.0.0

# Should output: v1.0.0
```

**If tag doesn't exist**:
```bash
# Create it
git tag -a v1.0.0 -m "Release v1.0.0 - Production Release"
git push origin v1.0.0
```

**Checkpoint**: Tag v1.0.0 exists on GitHub ✓

---

### ☐ Step 13: Create GitHub Release (3 minutes)

1. Go to: https://github.com/shandley/biometal/releases/new
2. Fill in:
   - **Choose a tag**: `v1.0.0` (select existing)
   - **Release title**: `v1.0.0 - Production Release`
   - **Description**: Copy from CHANGELOG.md (lines 8-149)
   - **Set as latest release**: ✓ (checked)
   - **Create a discussion**: ✓ (optional, but recommended)
3. Click "Publish release"

**What this does**: Triggers production workflow automatically!

**Checkpoint**: Release published, workflow starts ✓

---

### ☐ Step 14: Monitor Production Workflow (5-8 minutes)

1. Go to: https://github.com/shandley/biometal/actions
2. Click on "Publish to PyPI" workflow (should be running)
3. Watch jobs complete:
   - ⏳ Build wheels (4 platforms) - ~2-3 min each (parallel)
   - ⏳ Build source distribution - ~30s
   - ⏳ Publish to PyPI - ~1 min

**Expected**: All 6 jobs complete with green ✓

**If publish job waits**: Environment protection rules (approve in Actions)

**If any job fails**:
- Check error message
- Common issues:
  - Trusted publishing: Wait 2 min for propagation
  - Duplicate version: Cannot re-publish same version

**Checkpoint**: All jobs green, publish complete ✓

---

### ☐ Step 15: Verify on Production PyPI (2 minutes)

1. Go to: https://pypi.org/project/biometal-rs/
2. Check:
   - [ ] Package exists
   - [ ] Version is 1.0.0
   - [ ] Download files: 5 files
     - Mac ARM wheel
     - Mac x86_64 wheel
     - Linux x86_64 wheel
     - Linux ARM wheel
     - Source distribution (tar.gz)
   - [ ] Project description shows
   - [ ] Homepage link works
   - [ ] Documentation link works

**Checkpoint**: Package live on PyPI! ✓

---

### ☐ Step 16: Test Production Installation (2 minutes)

Open fresh terminal:

```bash
# Create test environment
python3 -m venv prod_test_env
source prod_test_env/bin/activate

# Install from PyPI (no --index-url needed!)
pip install biometal-rs==1.0.0

# Verify
python -c "import biometal; print(f'✓ Version: {biometal.__version__}')"

# Test operations
python << 'EOF'
import biometal
import gzip

# Create test data
with gzip.open('test.fq.gz', 'wt') as f:
    f.write('@read1\nATGCATGC\n+\nIIIIIIII\n')

# Test
stream = biometal.FastqStream.from_path('test.fq.gz')
record = next(stream)
gc = biometal.gc_content(bytes(record.sequence))
print(f'✓ GC content: {gc:.2%}')
print('✅ Production install works!')
EOF

# Cleanup
rm test.fq.gz
deactivate
rm -rf prod_test_env
```

**Expected**: Clean installation, all operations work

**Checkpoint**: Production pip install works! ✓

---

### 🎉 Phase 2 Complete!

**biometal v1.0.0 is now live on PyPI!**

Anyone in the world can now run:
```bash
pip install biometal-rs
```

---

## Phase 3: Post-Publication (5 minutes)

### ☐ Step 17: Announce Release

1. **GitHub Discussions**:
   - Go to: https://github.com/shandley/biometal/discussions
   - Create post: "biometal v1.0.0 Released!"
   - Content: Highlight features, link to docs

2. **Update README badges** (optional):
   - Add PyPI download count badge
   - Add version badge (auto-updates)

3. **Social media** (optional):
   - Twitter/X, LinkedIn, Reddit (r/bioinformatics)
   - "biometal v1.0.0: ARM-native bioinformatics, now on PyPI!"

---

### ☐ Step 18: Monitor for Issues

First 24 hours:
- [ ] Watch GitHub Issues
- [ ] Check PyPI download stats
- [ ] Monitor installation reports

Common first-day issues:
- Platform-specific installation problems
- Missing dependencies
- Documentation gaps

**Be ready to**: Hotfix if critical issue found (v1.0.1)

---

### ☐ Step 19: Update Documentation

If needed:
- [ ] Update QUICKSTART.md to emphasize `pip install`
- [ ] Add "Star on GitHub" badge
- [ ] Create CONTRIBUTORS.md if community contributions start

---

## ✅ Final Checklist

Mark when complete:

**TestPyPI (Validation)**:
- [ ] TestPyPI account created
- [ ] Trusted publishing configured (TestPyPI)
- [ ] GitHub `testpypi` environment created
- [ ] Test workflow ran successfully
- [ ] Package visible on test.pypi.org
- [ ] Test installation worked
- [ ] Test functionality passed

**Production PyPI**:
- [ ] PyPI account created (with 2FA)
- [ ] Trusted publishing configured (PyPI)
- [ ] GitHub `pypi` environment created
- [ ] Git tag v1.0.0 exists
- [ ] GitHub release published
- [ ] Production workflow completed
- [ ] Package live on pypi.org
- [ ] Production installation tested

**Post-Launch**:
- [ ] Release announced
- [ ] Monitoring in place
- [ ] Documentation updated

---

## 🆘 Troubleshooting

### Issue: "Trusted publishing not configured"

**Solution**:
1. Wait 2 minutes (propagation delay)
2. Verify configuration:
   - Workflow name matches exactly
   - Environment name matches exactly
   - Repository owner correct
3. Re-run workflow

### Issue: "File already exists"

**Cause**: Version already published to PyPI

**Solution**: Cannot delete from PyPI. Must bump version:
1. Update version in Cargo.toml and pyproject.toml
2. Commit and tag new version
3. Create new release

### Issue: Workflow fails on Linux ARM

**Cause**: Cross-compilation toolchain

**Status**: Already fixed in workflow (commit ac24b6f)

**If still fails**: Check logs for linker errors

### Issue: Can't trigger workflow

**Cause**: Not repository admin

**Solution**: Need "write" access to repository

---

## 📚 Resources

- **TestPyPI Guide**: [docs/TEST_PYPI_GUIDE.md](docs/TEST_PYPI_GUIDE.md)
- **PyPI Guide**: [docs/PYPI_PUBLISHING.md](docs/PYPI_PUBLISHING.md)
- **PyPI Help**: https://pypi.org/help/
- **Trusted Publishing**: https://docs.pypi.org/trusted-publishers/

---

## 🎯 Time Estimates

| Phase | Duration | Can Pause? |
|-------|----------|------------|
| TestPyPI Setup | 5 min | Yes |
| TestPyPI Test | 5 min | No (workflow running) |
| TestPyPI Verify | 5 min | Yes |
| **Total Test** | **15 min** | |
| Production Setup | 5 min | Yes |
| Production Publish | 5-8 min | No (workflow running) |
| Production Verify | 2 min | Yes |
| **Total Production** | **15 min** | |
| **GRAND TOTAL** | **30 minutes** | |

---

**Status**: Ready to begin
**Next Action**: Start with Step 1 (Create TestPyPI account)
**Goal**: `pip install biometal-rs` works worldwide! 🚀
