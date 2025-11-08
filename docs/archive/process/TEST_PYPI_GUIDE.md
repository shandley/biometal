# Testing PyPI Publishing with TestPyPI

Quick guide to safely test the PyPI publishing workflow using TestPyPI before production release.

---

## What is TestPyPI?

TestPyPI is a separate instance of the Python Package Index for testing purposes:
- URL: https://test.pypi.org
- Isolated from production PyPI
- Packages can be deleted
- Perfect for testing workflows

**Use it to validate**: Wheel builds, publishing process, installation

---

## Setup (5 minutes)

### 1. Create TestPyPI Account

1. Go to: https://test.pypi.org/account/register/
2. Create account (use same email as PyPI if you have one)
3. Verify email

### 2. Configure Trusted Publishing

1. Go to: https://test.pypi.org/manage/account/publishing/
2. Click "Add a new pending publisher"
3. Fill in:
   - **PyPI Project Name**: `biometal`
   - **Owner**: `shandley`
   - **Repository name**: `biometal`
   - **Workflow name**: `test-pypi-publish.yml`
   - **Environment name**: `testpypi`
4. Click "Add"

### 3. Create GitHub Environment

1. Go to: https://github.com/shandley/biometal/settings/environments
2. Click "New environment"
3. Name: `testpypi`
4. (Optional) Add protection rules
5. Click "Create environment"

---

## Running the Test

### Manual Trigger

1. Go to: https://github.com/shandley/biometal/actions/workflows/test-pypi-publish.yml
2. Click "Run workflow"
3. Select branch: `main`
4. Click "Run workflow"

The workflow will:
- Build wheels for Mac ARM and Linux x86_64 (faster test)
- Build source distribution
- Publish to TestPyPI

**Expected time**: 3-5 minutes

---

## What Gets Built

| Artifact | Size | Notes |
|----------|------|-------|
| Mac ARM wheel | ~1.0 MB | Python 3.11+ |
| Linux x86_64 wheel | ~0.9 MB | Python 3.11+ |
| Source dist | ~332 KB | All Python versions |

**Note**: Test workflow only builds 2 platforms (faster). Production builds 4.

---

## Monitoring the Build

### Watch Progress

1. Go to: https://github.com/shandley/biometal/actions
2. Click on the running workflow
3. Watch each job:
   - ✓ Build wheels (Mac ARM) - ~2 min
   - ✓ Build wheels (Linux x86_64) - ~2 min
   - ✓ Build sdist - ~30s
   - ✓ Test publish - ~30s

### Check Logs

Click on each job to see detailed logs:
- Build output
- Artifact uploads
- Publishing confirmation

---

## Verify on TestPyPI

### 1. Check Package Page

Go to: https://test.pypi.org/project/biometal-rs/

You should see:
- Version: 1.0.0
- Download files: 3 files (2 wheels + 1 sdist)
- Project description (from README.md)

### 2. Test Installation

```bash
# Install from TestPyPI
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple/ \
            biometal==1.0.0

# Note: --extra-index-url needed for dependencies
```

### 3. Verify Functionality

```python
import biometal

print(f"Version: {biometal.__version__}")

# Test basic functionality
stream = biometal.FastqStream.from_path("test.fq.gz")
for record in stream:
    gc = biometal.gc_content(bytes(record.sequence))
    print(f"{record.id}: GC = {gc:.2%}")
    break

print("✓ biometal works!")
```

---

## Troubleshooting

### Build Fails

**Check**:
- Rust toolchain installed on runner
- Python version correct (3.11)
- Cargo.toml and pyproject.toml versions match

**Look for**:
- Red ✗ in Actions tab
- Error messages in job logs

### Publishing Fails: "Trusted publishing not configured"

**Solution**:
1. Verify TestPyPI trusted publishing setup
2. Check workflow name: `test-pypi-publish.yml`
3. Check environment name: `testpypi`
4. Wait 1-2 minutes (propagation delay)

### Publishing Fails: "File already exists"

**Not a problem**: TestPyPI uses `skip-existing: true`

**Explanation**: Re-running workflow with same version is safe

### Can't Install: "No matching distribution"

**Solution**: Add `--extra-index-url https://pypi.org/simple/`

**Reason**: Dependencies aren't on TestPyPI, need main PyPI

---

## Success Criteria

✅ **All checks passed**:
- [ ] Workflow completed without errors
- [ ] All 3 artifacts built (2 wheels + sdist)
- [ ] Package visible on test.pypi.org
- [ ] `pip install` works from TestPyPI
- [ ] `import biometal` works
- [ ] Basic operations work (streaming, GC content)

**If all ✓**: Ready for production PyPI!

---

## Cleanup (Optional)

TestPyPI packages can be deleted:

1. Go to: https://test.pypi.org/manage/project/biometal-rs/releases/
2. Select version
3. Click "Options" → "Delete"

**Note**: Not necessary, TestPyPI is for testing

---

## Next Steps After Successful Test

1. **Configure Production PyPI** (if not already done):
   - URL: https://pypi.org/manage/account/publishing/
   - Setup same as TestPyPI but use `publish-pypi.yml`

2. **Create v1.0.0 Release**:
   - Tag: `v1.0.0`
   - Triggers production workflow
   - Publishes to pypi.org

3. **Verify Production**:
   ```bash
   pip install biometal-rs==1.0.0  # No --index-url needed
   ```

---

## Differences: Test vs Production

| Feature | Test Workflow | Production Workflow |
|---------|--------------|---------------------|
| Trigger | Manual only | GitHub Release |
| Platforms | 2 (faster) | 4 (complete) |
| Target | TestPyPI | PyPI |
| URL | test.pypi.org | pypi.org |
| Skip existing | Yes | No (fail if exists) |

---

## Tips

### 1. Test Early, Test Often

Run test workflow on feature branches before merging

### 2. Verify Platform Coverage

Production builds 4 platforms. Test builds 2 (representative sample).

### 3. Check Installation

Always test `pip install` from TestPyPI before production

### 4. Use Test for Breaking Changes

Major version bumps? Test on TestPyPI first.

---

## Example Test Session

```bash
# 1. Trigger workflow
# Go to GitHub Actions → test-pypi-publish.yml → Run workflow

# 2. Wait for completion (~5 minutes)

# 3. Install from TestPyPI
python3 -m pip install --index-url https://test.pypi.org/simple/ \
                       --extra-index-url https://pypi.org/simple/ \
                       biometal==1.0.0

# 4. Test functionality
python3 -c "
import biometal
print(f'✓ Version: {biometal.__version__}')

stream = biometal.FastqStream.from_path('test.fq.gz')
record = next(stream)
gc = biometal.gc_content(bytes(record.sequence))
print(f'✓ GC content: {gc:.2%}')
print('✓ All tests passed!')
"

# 5. Uninstall (optional)
pip uninstall biometal -y
```

---

## FAQ

### Q: Do I need to test every time?

**A**: No. Test once to validate workflow. Re-test for:
- Major version changes
- Workflow modifications
- New platforms added

### Q: Can I delete TestPyPI packages?

**A**: Yes. Unlike PyPI, TestPyPI allows deletion.

### Q: How long do TestPyPI packages stay?

**A**: No expiration, but may be purged during maintenance.

### Q: Do downloads from TestPyPI count?

**A**: No. TestPyPI downloads aren't tracked in PyPI statistics.

---

## Resources

- **TestPyPI**: https://test.pypi.org
- **TestPyPI Help**: https://test.pypi.org/help/
- **Trusted Publishing**: https://docs.pypi.org/trusted-publishers/

---

**Status**: Ready to test
**Next**: Run test workflow → Verify → Publish to production
