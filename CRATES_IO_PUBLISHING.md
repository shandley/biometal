# Publishing biometal v1.0.0 to crates.io

**Status**: ✅ Ready to publish!
**Estimated Time**: 3 minutes

---

## ✅ Pre-Publishing Complete

I've already completed:
- ✅ Updated author email to `handley.scott@gmail.com`
- ✅ Ran `cargo publish --dry-run` - **SUCCESS**
- ✅ Verified package contents (84 files, 812.1 KB)
- ✅ Committed and pushed changes
- ✅ All 121 tests passing

**The package is ready to publish!**

---

## 🚀 Two Steps to Complete

### Step 1: Get Your crates.io API Token (1 minute)

1. **Go to**: https://crates.io/settings/tokens
2. **Log in** with your GitHub account (if not already logged in)
3. **Click**: "New Token"
4. **Name**: `biometal-publishing` (or any name you like)
5. **Permissions**: Leave default (publish new crates)
6. **Click**: "Generate"
7. **Copy** the token (starts with `cio_...`)

**Important**: Save this token somewhere safe! You won't be able to see it again.

---

### Step 2: Publish to crates.io (2 minutes)

Open your terminal and run these commands:

```bash
# Login to crates.io (one-time setup)
cargo login <paste-your-token-here>

# Publish biometal v1.0.0
cargo publish
```

**Expected output**:
```
   Updating crates.io index
   Packaging biometal v1.0.0
   Packaged 84 files, 812.1KiB (334.4KiB compressed)
   Verifying biometal v1.0.0
   Compiling biometal v1.0.0
    Finished `dev` profile target(s) in 8.62s
   Uploading biometal v1.0.0
✅ Successfully uploaded biometal v1.0.0 to crates.io
```

**Time**: ~10 seconds for upload, ~5 minutes for indexing

---

## 🎉 Verification

Once published, verify at:
- **Crate page**: https://crates.io/crates/biometal
- **Documentation**: https://docs.rs/biometal (auto-generated in ~10 minutes)

Test installation:
```bash
cargo new test_biometal
cd test_biometal
cargo add biometal
cargo build
```

---

## 📦 What Gets Published

Your package includes:
- **84 files** (812.1 KB total, 334.4 KB compressed)
- Source code (src/, benches/, examples/)
- Documentation (README.md, CLAUDE.md, docs/)
- Licenses (LICENSE-APACHE, LICENSE-MIT)
- Metadata (Cargo.toml, CHANGELOG.md)

**Network feature**: Enabled by default (as configured)

---

## 🔄 Future Updates

For future versions (v1.0.1, etc.):

1. Update version in `Cargo.toml`:
   ```toml
   version = "1.0.1"
   ```

2. Update `CHANGELOG.md` with changes

3. Commit and tag:
   ```bash
   git commit -am "chore: Bump version to 1.0.1"
   git tag v1.0.1
   git push && git push --tags
   ```

4. Publish:
   ```bash
   cargo publish
   ```

**No need to re-login** - your token is saved locally.

---

## 🆘 Troubleshooting

### Issue: "error: failed to authenticate to registry"

**Solution**: Run `cargo login <your-token>` again

---

### Issue: "crate version `1.0.0` is already uploaded"

**Cause**: Version already published (can't overwrite)

**Solution**: Bump version to 1.0.1 and republish

---

### Issue: "error: 1 files in the working directory contain changes"

**Solution**: Commit changes first, or use `--allow-dirty` (not recommended)

---

### Issue: "error: no upload token found"

**Solution**: You need to run `cargo login` first

---

## 📊 Post-Publishing

### Update README Badges (Optional)

Add to the top of README.md:
```markdown
[![Crates.io](https://img.shields.io/crates/v/biometal.svg)](https://crates.io/crates/biometal)
[![Documentation](https://docs.rs/biometal/badge.svg)](https://docs.rs/biometal)
[![Downloads](https://img.shields.io/crates/d/biometal.svg)](https://crates.io/crates/biometal)
```

### Monitor Stats

- **Downloads**: https://crates.io/crates/biometal/stats
- **Reverse dependencies**: See who's using your crate
- **docs.rs build**: Monitor at https://docs.rs/releases/queue

### Announce

**Rust community**:
- Reddit r/rust - "Show & Tell" thread
- Twitter/X with #rustlang hashtag
- This Week in Rust newsletter submission

**Example announcement**:
> 🎉 biometal v1.0.0 is now on crates.io!
>
> ARM-native bioinformatics with 16-25× NEON speedup. Analyze 5TB datasets on consumer laptops with constant ~5 MB memory.
>
> Add to your project: `cargo add biometal`
>
> https://crates.io/crates/biometal

---

## 🎯 Success Checklist

Mark off when complete:

- [ ] Got crates.io API token
- [ ] Ran `cargo login`
- [ ] Ran `cargo publish`
- [ ] Verified at https://crates.io/crates/biometal
- [ ] Tested `cargo add biometal` in new project
- [ ] (Optional) Updated README badges
- [ ] (Optional) Announced to Rust community

---

## 📚 Resources

- **crates.io Publishing Guide**: https://doc.rust-lang.org/cargo/reference/publishing.html
- **API Tokens**: https://crates.io/settings/tokens
- **docs.rs**: https://docs.rs/about
- **Cargo Book**: https://doc.rust-lang.org/cargo/

---

## 🎊 What This Achieves

Once published, Rust developers worldwide can:

```toml
[dependencies]
biometal = "1.0"
```

Instead of:
```toml
[dependencies]
biometal = { git = "https://github.com/shandley/biometal" }
```

**Better experience**: Semantic versioning, easy updates, discoverable via search.

---

**Ready?** Start with Step 1 (Get API Token) above!

**Questions?** Check troubleshooting or reach out via GitHub Issues.

---

**Last Updated**: November 5, 2025
**Version**: 1.0.0
**Package Size**: 334.4 KB (compressed)
