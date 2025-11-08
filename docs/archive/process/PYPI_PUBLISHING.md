# Publishing biometal to PyPI

Guide for publishing biometal Python bindings to the Python Package Index (PyPI).

---

## Prerequisites

### 1. PyPI Account Setup

**Create accounts** (if you don't have them):
- Production: https://pypi.org/account/register/
- Test (optional): https://test.pypi.org/account/register/

**Enable 2FA**: PyPI requires two-factor authentication for package maintainers.

### 2. GitHub Setup

**Configure trusted publishing** (recommended, no API tokens needed):

1. Go to https://pypi.org/manage/account/publishing/
2. Add a new publisher:
   - **PyPI Project Name**: `biometal`
   - **Owner**: `shandley`
   - **Repository name**: `biometal`
   - **Workflow name**: `publish-pypi.yml`
   - **Environment name**: `pypi`

3. Create PyPI environment in GitHub:
   - Go to: https://github.com/shandley/biometal/settings/environments
   - Click "New environment"
   - Name: `pypi`
   - Add protection rules (optional): Require review before deployment

### 3. Local Tools

Install required tools:

```bash
# Install maturin (Rust → Python bridge)
pip install maturin

# Install twine (for manual uploads, optional)
pip install twine

# Verify Rust toolchain
rustc --version  # Need 1.70+
```

---

## Quick Publish (Automated)

### Via GitHub Release

1. **Tag the release**:
   ```bash
   git tag v1.0.0
   git push origin v1.0.0
   ```

2. **Create GitHub release**:
   - Go to: https://github.com/shandley/biometal/releases/new
   - Select tag: `v1.0.0`
   - Release title: `v1.0.0 - Production Release`
   - Description: Copy from CHANGELOG.md
   - Click "Publish release"

3. **Automated build & publish**:
   - GitHub Actions triggers `.github/workflows/publish-pypi.yml`
   - Builds wheels for: macOS ARM, macOS x86_64, Linux x86_64, Linux ARM
   - Publishes to PyPI via trusted publishing

4. **Verify**:
   - Check workflow: https://github.com/shandley/biometal/actions
   - Check PyPI: https://pypi.org/project/biometal-rs/

**Done!** Users can now `pip install biometal-rs`.

---

## Manual Publish (Local Build)

### 1. Build Wheels

```bash
# Build for current platform
maturin build --release --features python

# Output: dist/biometal-1.0.0-*.whl
```

**Cross-platform builds**:

```bash
# macOS ARM (if on macOS ARM)
maturin build --release --features python --target aarch64-apple-darwin

# macOS x86_64 (if on macOS Intel)
maturin build --release --features python --target x86_64-apple-darwin

# Linux x86_64 (if on Linux)
maturin build --release --features python --target x86_64-unknown-linux-gnu
```

### 2. Build Source Distribution

```bash
maturin sdist

# Output: dist/biometal-1.0.0.tar.gz
```

### 3. Test Locally

```bash
# Install from wheel
pip install dist/biometal-1.0.0-*.whl

# Test
python -c "import biometal; print(biometal.__version__)"

# Run test suite
python test_python_bindings.py
```

### 4. Upload to PyPI

**Option 1: Using maturin (recommended)**:
```bash
maturin publish --username __token__ --password $PYPI_TOKEN
```

**Option 2: Using twine**:
```bash
twine upload dist/*
```

**For test PyPI** (practice first):
```bash
maturin publish --repository testpypi --username __token__ --password $TESTPYPI_TOKEN
```

---

## Version Management

### Before Each Release

1. **Update versions**:
   ```bash
   # Update Cargo.toml
   sed -i '' 's/version = "0.2.3"/version = "1.0.0"/' Cargo.toml

   # Update pyproject.toml
   sed -i '' 's/version = "0.2.3"/version = "1.0.0"/' pyproject.toml

   # Verify
   grep version Cargo.toml pyproject.toml
   ```

2. **Update CHANGELOG.md**:
   - Document all changes
   - Link to evidence (experiments)
   - Note breaking changes

3. **Update README.md**:
   - Verify installation instructions
   - Update status section

4. **Run tests**:
   ```bash
   # Rust tests
   cargo test --all-features

   # Python tests
   maturin develop --release --features python
   python test_python_bindings.py

   # Clippy
   cargo clippy --all-features -- -D warnings
   ```

5. **Build and test locally**:
   ```bash
   maturin build --release --features python
   pip install --force-reinstall dist/biometal-1.0.0-*.whl
   python test_python_bindings.py
   ```

6. **Commit and tag**:
   ```bash
   git add Cargo.toml pyproject.toml CHANGELOG.md README.md
   git commit -m "chore: Bump version to 1.0.0"
   git tag v1.0.0
   git push && git push --tags
   ```

---

## Troubleshooting

### Build Fails: "cargo not found"

**Solution**: Install Rust toolchain:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

---

### Build Fails: "linker 'cc' not found"

**Solution**: Install C compiler:

**macOS**:
```bash
xcode-select --install
```

**Ubuntu**:
```bash
sudo apt install build-essential
```

---

### Upload Fails: "Invalid or non-existent authentication"

**Solution**: Use API token, not password:

1. Generate token: https://pypi.org/manage/account/token/
2. Scope: "Entire account" or "Project: biometal"
3. Copy token (starts with `pypi-`)
4. Use in upload:
   ```bash
   maturin publish --username __token__ --password pypi-AgENd...
   ```

---

### Upload Fails: "Filename already exists"

**Solution**: Can't overwrite existing version. Bump version:
```bash
# Update version in Cargo.toml and pyproject.toml
# Then rebuild and upload
```

---

### Wheel Fails to Install: "Unsupported platform"

**Solution**: Build for correct platform. Common platforms:

| Platform | Target |
|----------|--------|
| Mac ARM | `aarch64-apple-darwin` |
| Mac Intel | `x86_64-apple-darwin` |
| Linux x86_64 | `x86_64-unknown-linux-gnu` |
| Linux ARM | `aarch64-unknown-linux-gnu` |

```bash
rustup target add aarch64-apple-darwin
maturin build --release --features python --target aarch64-apple-darwin
```

---

## Platform Support Strategy

### Primary Targets (Must Build)

1. **macOS ARM** (Mac M1/M2/M3/M4)
   - Target: `aarch64-apple-darwin`
   - Optimized with NEON (16-25× speedup)
   - Primary democratization platform

2. **Linux x86_64** (Most servers)
   - Target: `x86_64-unknown-linux-gnu`
   - Scalar fallback (portable)

### Secondary Targets (Nice to Have)

3. **macOS x86_64** (Intel Macs)
   - Target: `x86_64-apple-darwin`
   - Scalar fallback

4. **Linux ARM** (Graviton, Raspberry Pi)
   - Target: `aarch64-unknown-linux-gnu`
   - NEON support (6-10× speedup)

### Not Currently Supported

- **Windows**: Not tested, community contributions welcome
- **32-bit platforms**: Not supported

---

## Release Checklist

Use this checklist for each release:

- [ ] **Update versions** (Cargo.toml, pyproject.toml)
- [ ] **Update CHANGELOG.md** (features, fixes, breaking changes)
- [ ] **Run all tests** (`cargo test`, `python test_python_bindings.py`)
- [ ] **Run clippy** (`cargo clippy --all-features`)
- [ ] **Build locally** (`maturin build --release --features python`)
- [ ] **Test local wheel** (install and run tests)
- [ ] **Commit changes** (version bump, changelog)
- [ ] **Create git tag** (`git tag v1.0.0`)
- [ ] **Push with tags** (`git push && git push --tags`)
- [ ] **Create GitHub release** (triggers automated build)
- [ ] **Verify build** (check GitHub Actions)
- [ ] **Verify PyPI** (check package appears)
- [ ] **Test pip install** (`pip install biometal-rs==1.0.0`)
- [ ] **Announce** (GitHub discussions, Twitter, etc.)

---

## Testing Before Release

### Test on TestPyPI First

1. **Build and upload to TestPyPI**:
   ```bash
   maturin publish --repository testpypi
   ```

2. **Install from TestPyPI**:
   ```bash
   pip install --index-url https://test.pypi.org/simple/ biometal==1.0.0
   ```

3. **Test functionality**:
   ```python
   import biometal
   print(biometal.__version__)

   stream = biometal.FastqStream.from_path("test.fq.gz")
   for record in stream:
       print(record.id)
       break
   ```

4. **If successful**, proceed to production PyPI.

---

## Post-Release Tasks

### 1. Update Documentation

- [ ] Verify README.md installation instructions work
- [ ] Update [QUICKSTART.md](../QUICKSTART.md) if needed
- [ ] Update [PYTHON.md](../PYTHON.md) if API changed

### 2. Announce Release

- [ ] GitHub Discussions post
- [ ] Twitter/X announcement
- [ ] Relevant bioinformatics forums
- [ ] Email collaborators/early adopters

### 3. Monitor Issues

- [ ] Watch for installation issues
- [ ] Watch for platform-specific problems
- [ ] Respond to bug reports quickly

### 4. Update Conda (Future)

Once on PyPI, create conda-forge recipe:
- Submit PR to https://github.com/conda-forge/staged-recipes
- Follow conda-forge guidelines

---

## Continuous Integration

### Automated Testing (on PR)

The `.github/workflows/ci.yml` should test:
- [ ] Rust tests on multiple platforms
- [ ] Python tests on multiple platforms
- [ ] Clippy linting
- [ ] Format checking (`cargo fmt --check`)

### Automated Publishing (on Release)

The `.github/workflows/publish-pypi.yml` handles:
- [ ] Building wheels for all platforms
- [ ] Building source distribution
- [ ] Publishing to PyPI via trusted publishing

---

## Security Best Practices

### 1. Use Trusted Publishing

**Recommended**: Use GitHub Actions with trusted publishing (no API tokens in repo).

**Setup**: Configure at https://pypi.org/manage/account/publishing/

### 2. If Using API Tokens

**DO**:
- Use GitHub Secrets for tokens
- Scope tokens to single project
- Rotate tokens regularly

**DON'T**:
- Commit tokens to repo
- Share tokens in public
- Use account-wide tokens for automation

### 3. Code Signing

**Future**: Sign releases with GPG:
```bash
git tag -s v1.0.0 -m "Release 1.0.0"
```

---

## Resources

- **PyPI**: https://pypi.org/project/biometal-rs/
- **Maturin Docs**: https://www.maturin.rs/
- **PyPI Publishing Guide**: https://packaging.python.org/tutorials/packaging-projects/
- **Trusted Publishing**: https://docs.pypi.org/trusted-publishers/

---

## Contact

**Questions or issues**:
- GitHub Issues: https://github.com/shandley/biometal/issues
- GitHub Discussions: https://github.com/shandley/biometal/discussions

---

**Last Updated**: November 5, 2025 (v1.0.0)
