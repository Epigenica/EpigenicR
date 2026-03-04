# Git Setup Guide for EpigenicR

## Prerequisites

1. GitHub account with access to Epigenica organization
2. Git configured locally with your credentials

## Step-by-Step Guide

### 1. Create Remote Repository on GitHub

**Option A: Via GitHub Web Interface (Recommended for Organizations)**

1. Go to https://github.com/epigenica
2. Click "New repository" (or go to https://github.com/organizations/epigenica/repositories/new)
3. Repository settings:
   - **Name**: `EpigenicR`
   - **Description**: `R package for EpiFinder platform data analysis - Tools for epigenomic data analysis and visualization`
   - **Visibility**: ✅ Private (for internal use)
   - **Initialize**: ❌ DO NOT add README, .gitignore, or license (we already have these)
4. Click "Create repository"
5. Copy the repository URL: `git@github.com:epigenica/EpigenicR.git` (SSH) or `https://github.com/epigenica/EpigenicR.git` (HTTPS)

**Option B: Via GitHub CLI (if installed)**

```bash
gh repo create epigenica/EpigenicR --private --description "R package for EpiFinder platform data analysis" --source=. --remote=origin
```

### 2. Update .gitignore (Important!)

The current .gitignore is minimal. Run this to add standard R package ignores:

```bash
cat >> .gitignore << 'EOF'

# R package build artifacts
*.tar.gz
*.tgz
/check/
/revdep/
/*.Rcheck/

# Documentation
/doc/
/Meta/

# pkgdown site
docs/
_pkgdown.yml

# Session Data
.RData
.Rhistory

# Temporary files
*~
*.swp
*.swo

# OS files
.DS_Store
Thumbs.db

# IDE
.vscode/
.idea/

# Sensitive files (if any)
.Renviron
.env
EOF
```

### 3. Initialize and Commit Everything

```bash
# Check what files will be committed
git status

# Add all files
git add .

# Review what's staged (optional but recommended)
git status

# Make initial commit
git commit -m "Initial commit: EpigenicR package v0.1.0

- Add core functions for epigenomic data analysis
- Include toy dataset for testing
- Add comprehensive documentation and README
- Set up GPL-3 license
- Include example scripts and QC plotting functions"

# Verify commit
git log --oneline
```

### 4. Add Remote Repository

Replace `<REPO_URL>` with your actual repository URL from Step 1:

```bash
# SSH (recommended if you have SSH keys set up)
git remote add origin git@github.com:epigenica/EpigenicR.git

# OR HTTPS (if you prefer token authentication)
# git remote add origin https://github.com/epigenica/EpigenicR.git

# Verify remote was added
git remote -v
```

### 5. Push to GitHub

```bash
# Push to remote (first time)
git push -u origin master

# Or if using 'main' as default branch:
# git branch -M main
# git push -u origin main
```

### 6. Verify on GitHub

1. Go to https://github.com/epigenica/EpigenicR
2. Verify all files are there
3. Check that README renders correctly with images
4. Verify repository is marked as "Private"

## Post-Setup: Team Access

### Set Team Permissions

1. Go to repository Settings → Manage access
2. Add team members or teams:
   - **Admin**: Development leads
   - **Write**: Active developers
   - **Read**: Other team members

### Create Initial Release (Optional)

```bash
# Tag the initial version
git tag -a v0.1.0 -m "Initial release: EpigenicR v0.1.0"
git push origin v0.1.0
```

Then create a GitHub release:
1. Go to repository → Releases
2. Click "Draft a new release"
3. Choose tag: v0.1.0
4. Title: "EpigenicR v0.1.0 - Initial Release"
5. Add release notes
6. Publish release

## Installation Instructions for Team

After pushing, team members can install with:

```r
# Using HTTPS with personal access token
remotes::install_github("epigenica/EpigenicR", 
                        auth_token = "your_github_token")

# Or using SSH (if SSH keys configured)
remotes::install_github("epigenica/EpigenicR")
```

## Troubleshooting

### Authentication Issues

**For HTTPS:**
- Create a Personal Access Token: GitHub Settings → Developer settings → Personal access tokens
- Use token instead of password when prompted

**For SSH:**
- Generate SSH key: `ssh-keygen -t ed25519 -C "your_email@epigenica.com"`
- Add to GitHub: GitHub Settings → SSH and GPG keys

### Large Files Warning

If you get warnings about large files (>50MB):
- BigWig files in `inst/extdata/toy_dataset/` are ~13MB total (should be fine)
- If needed, use Git LFS for large files:
  ```bash
  git lfs install
  git lfs track "*.bw"
  git add .gitattributes
  ```

### Branch Name Conflicts

If GitHub default is 'main' but local is 'master':
```bash
git branch -M main
git push -u origin main
```

## Quick Reference Commands

```bash
# Check status
git status

# View remote
git remote -v

# Pull latest changes
git pull origin master

# Push changes
git push origin master

# Create new branch
git checkout -b feature/new-feature

# Switch branches
git checkout master
```

## Next Steps

1. ✅ Push code to GitHub
2. Set up branch protection rules (Settings → Branches)
3. Add CODEOWNERS file for automatic review requests
4. Set up GitHub Actions for R CMD check (optional)
5. Create project board for issue tracking
