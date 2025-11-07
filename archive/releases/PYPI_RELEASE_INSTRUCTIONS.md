# PyPI Production Release - Manual Steps

**Status**: Ready to publish biometal-rs v1.0.0 to PyPI
**Estimated Time**: 15-20 minutes

---

## ✅ Pre-Release Verification (COMPLETED)

- ✅ Git tag v1.0.0 exists
- ✅ Versions match (Cargo.toml: 1.0.0, pyproject.toml: 1.0.0)
- ✅ README references correct package name (biometal-rs)
- ✅ TestPyPI publishing successful
- ✅ Release notes prepared
- ✅ All 121 tests passing
- ✅ GitHub workflows configured

---

## 🔧 Manual Steps Required

### Step 1: Create Production PyPI Account (5 minutes)

**If you don't already have one:**

1. Go to: https://pypi.org/account/register/
2. Fill in:
   - Username: [your choice]
   - Email: [your email]
   - Password: [secure password]
3. Verify email (check inbox)
4. **REQUIRED**: Enable 2FA
   - Go to: https://pypi.org/manage/account/
   - Click "Add 2FA with authentication application"
   - Scan QR code with authenticator app (Google Authenticator, Authy, etc.)
   - Enter verification code

**Checkpoint**: You can log in to pypi.org with 2FA ✓

---

### Step 2: Configure PyPI Trusted Publishing (5 minutes)

1. Log in to PyPI: https://pypi.org
2. Go to: https://pypi.org/manage/account/publishing/
3. Click "Add a new pending publisher"
4. Fill in the form **EXACTLY** as shown:

   ```
   PyPI Project Name: biometal-rs
   Owner: shandley
   Repository name: biometal
   Workflow name: publish-pypi.yml
   Environment name: pypi
   ```

5. Click "Add"

**Important Notes**:
- Project name is `biometal-rs` (not `biometal`)
- Workflow name is `publish-pypi.yml` (not `test-pypi-publish.yml`)
- Environment name is `pypi` (not `testpypi`)

**Checkpoint**: You see "Pending publisher added for biometal-rs" ✓

---

### Step 3: Create GitHub `pypi` Environment (3 minutes)

1. Go to: https://github.com/shandley/biometal/settings/environments
2. Click "New environment" button
3. Name: `pypi` (exactly, case-sensitive)
4. **Recommended**: Add protection rules for safety:
   - ✓ "Required reviewers" → Add yourself
   - This gives you a final approval gate before publishing
5. Click "Configure environment"

**Checkpoint**: You see `pypi` in environments list ✓

---

### Step 4: Create GitHub Release (2 minutes)

1. Go to: https://github.com/shandley/biometal/releases/new

2. Fill in the form:
   - **Choose a tag**: Select `v1.0.0` from dropdown (already exists)
   - **Release title**: `v1.0.0 - Production Release`
   - **Description**: Copy contents from `RELEASE_NOTES_v1.0.0.md`
   - **Set as the latest release**: ✓ (checked)
   - **Create a discussion for this release**: ✓ (optional, recommended)

3. Click "Publish release" button

**What happens**:
- GitHub release is created
- Workflow `.github/workflows/publish-pypi.yml` triggers automatically
- Builds wheels for 4 platforms
- Publishes to PyPI

**Checkpoint**: Release published, workflow starts running ✓

---

### Step 5: Monitor Workflow (5-8 minutes)

1. Go to: https://github.com/shandley/biometal/actions
2. Click on "Publish to PyPI" workflow (should show as running)
3. Watch jobs complete:

   ```
   ⏳ Build wheels (macos-aarch64)    - ~2-3 min
   ⏳ Build wheels (macos-x86_64)     - ~2-3 min
   ⏳ Build wheels (linux-x86_64)     - ~2-3 min
   ⏳ Build wheels (linux-aarch64)    - ~2-3 min
   ⏳ Build source distribution        - ~30s
   ⏳ Publish to PyPI                  - ~1 min
   ```

**If "Publish to PyPI" job waits for approval**:
- This is normal if you set up environment protection
- Click "Review deployments" → Approve

**Expected Result**: All 6 jobs complete with green checkmarks ✓

---

### Step 6: Verify on PyPI (2 minutes)

1. Go to: https://pypi.org/project/biometal-rs/
2. Verify:
   - [ ] Package exists
   - [ ] Version shows as 1.0.0
   - [ ] "Download files" shows 5 files:
     - `biometal_rs-1.0.0-cp311-cp311-macosx_11_0_arm64.whl`
     - `biometal_rs-1.0.0-cp311-cp311-macosx_10_12_x86_64.whl`
     - `biometal_rs-1.0.0-cp311-cp311-manylinux_2_17_x86_64.whl`
     - `biometal_rs-1.0.0-cp311-cp311-manylinux_2_17_aarch64.whl`
     - `biometal-rs-1.0.0.tar.gz`
   - [ ] Project description renders correctly (from README.md)
   - [ ] Links work (Homepage, Documentation, Repository)

**Checkpoint**: Package is live on PyPI! ✓

---

### Step 7: Test Installation (3 minutes)

Open a terminal and run:

```bash
# Create fresh test environment
python3 -m venv pypi_test_env
source pypi_test_env/bin/activate

# Install from PyPI (production!)
pip install biometal-rs==1.0.0

# Verify version
python -c "import biometal; print(f'✓ Version: {biometal.__version__}')"
```

**Expected output**: `✓ Version: 1.0.0`

**Test basic functionality**:

```bash
python << 'EOF'
import biometal
import gzip

# Create test FASTQ
with gzip.open('test.fq.gz', 'wt') as f:
    f.write('@read1\nATGCATGC\n+\nIIIIIIII\n')

# Test streaming
stream = biometal.FastqStream.from_path('test.fq.gz')
record = next(stream)

# Test operations
gc = biometal.gc_content(bytes(record.sequence))
counts = biometal.count_bases(bytes(record.sequence))
mean_q = biometal.mean_quality(bytes(record.quality))

print(f'✓ Streaming: {record.id}')
print(f'✓ GC content: {gc:.2%}')
print(f'✓ Base counts: {counts}')
print(f'✓ Mean quality: {mean_q:.1f}')
print('\n✅ ALL TESTS PASSED - biometal-rs is live!')
EOF
```

**Expected output**:
```
✓ Streaming: read1
✓ GC content: 50.00%
✓ Base counts: {'A': 2, 'C': 2, 'G': 2, 'T': 2}
✓ Mean quality: 40.0

✅ ALL TESTS PASSED - biometal-rs is live!
```

**Cleanup**:
```bash
rm test.fq.gz
deactivate
rm -rf pypi_test_env
```

**Checkpoint**: Production installation works perfectly! ✓

---

## 🎉 Success!

**biometal-rs v1.0.0 is now live on PyPI!**

Anyone in the world can now run:
```bash
pip install biometal-rs
```

---

## 📣 Post-Release Actions (Optional)

### 1. Announce the Release

**GitHub Discussions** (recommended):
- Go to: https://github.com/shandley/biometal/discussions
- Category: Announcements
- Title: "biometal v1.0.0 Released - Now on PyPI!"
- Content:
  ```markdown
  🎉 biometal v1.0.0 is now available on PyPI!

  Install with:
  ```bash
  pip install biometal-rs
  ```

  Key features:
  - 16-25× ARM NEON speedup
  - Streaming architecture (constant ~5 MB memory)
  - Network streaming (analyze without downloading)
  - Python 3.9-3.14 support

  Full release notes: https://github.com/shandley/biometal/releases/tag/v1.0.0
  ```

**Social Media** (if desired):
- Twitter/X, LinkedIn, Mastodon
- Reddit: r/bioinformatics, r/Python
- Message: "biometal v1.0.0: ARM-native bioinformatics library now on PyPI. Analyze 5TB datasets on laptops with 16-25× NEON speedup. #bioinformatics #python #opensource"

### 2. Update Badges (optional)

Add to README.md:
```markdown
[![PyPI version](https://badge.fury.io/py/biometal-rs.svg)](https://pypi.org/project/biometal-rs/)
[![Downloads](https://pepy.tech/badge/biometal-rs)](https://pepy.tech/project/biometal-rs)
```

### 3. Monitor for Issues

**First 24-48 hours**:
- [ ] Watch GitHub Issues for installation problems
- [ ] Check PyPI download stats: https://pypistats.org/packages/biometal-rs
- [ ] Monitor for platform-specific issues

**Common first-day issues**:
- Platform compatibility (usually x86_64 vs ARM)
- Python version mismatches
- Missing system dependencies

---

## 🆘 Troubleshooting

### Issue: "Trusted publishing not configured"

**Symptoms**: Workflow fails at "Publish to PyPI" step

**Solution**:
1. Verify configuration at https://pypi.org/manage/account/publishing/
2. Check all fields match exactly (especially `biometal-rs` not `biometal`)
3. Wait 2 minutes for propagation
4. Re-run workflow from GitHub Actions UI

---

### Issue: "File already exists"

**Symptoms**: Upload fails with "filename has already been used"

**Cause**: Version 1.0.0 already published (cannot overwrite on PyPI)

**Solution**:
- This is permanent - PyPI doesn't allow re-uploading same version
- If needed, bump to 1.0.1 and release again
- Prevention: Always test on TestPyPI first

---

### Issue: Workflow builds but doesn't publish

**Symptoms**: Build jobs succeed, publish job doesn't start

**Cause**: Environment protection rules waiting for approval

**Solution**:
1. Go to workflow run page
2. Look for yellow "Review deployments" button
3. Click it and approve

---

### Issue: Can't install on specific platform

**Symptoms**: `pip install biometal-rs` fails with "no matching distribution"

**Cause**: Platform not supported or wheel didn't build

**Solution**:
1. Check which wheels exist on PyPI
2. Verify platform: `python -c "import platform; print(platform.machine())"`
3. Fall back to source install: `pip install biometal-rs --no-binary :all:`

---

## 📚 Reference Documentation

- **PYPI_SETUP_CHECKLIST.md** - Complete 30-minute guide
- **docs/PYPI_PUBLISHING.md** - Detailed publishing reference
- **docs/TEST_PYPI_GUIDE.md** - TestPyPI testing guide
- **CHANGELOG.md** - Version history

---

## ✅ Final Checklist

Mark off as you complete:

**Pre-Release** (Automated - DONE):
- [x] Git tag v1.0.0 exists
- [x] Versions match (1.0.0)
- [x] README references biometal-rs
- [x] TestPyPI successful
- [x] Release notes prepared

**Production Publishing** (Manual - YOUR STEPS):
- [ ] PyPI account created with 2FA
- [ ] Trusted publishing configured (biometal-rs)
- [ ] GitHub `pypi` environment created
- [ ] GitHub release published
- [ ] Workflow completed successfully
- [ ] Package visible on pypi.org
- [ ] `pip install biometal-rs` tested
- [ ] Basic functionality verified

**Post-Release** (Optional):
- [ ] Release announced
- [ ] Monitoring in place
- [ ] Badges updated

---

**Ready to Begin?**

Start with **Step 1** (Create PyPI Account) above.

**Estimated Total Time**: 15-20 minutes

**Questions?** Check troubleshooting section or the full guides in `docs/`

---

**Good luck! 🚀**
