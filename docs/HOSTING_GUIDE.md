# HPC-ScaleTest Documentation - Hosting Guide

## 📦 What You Have

I've created a complete Sphinx documentation package for your HPC-ScaleTest project that includes:

- ✅ Full Sphinx setup with `conf.py` configuration
- ✅ Comprehensive documentation pages (getting started, installation, user guide, etc.)
- ✅ API reference structure
- ✅ Examples and tutorials
- ✅ Build automation (Makefile)
- ✅ Ready-to-deploy configuration files

## 🚀 Where to Host Your Documentation

### Option 1: Read the Docs (RECOMMENDED) ⭐

**Best for:** Open source projects, automatic builds, free hosting

**Why Choose This:**
- 🆓 Free for open source projects
- 🔄 Automatic builds on every commit
- 📚 Version management (latest, stable, v1.0.0, etc.)
- 🔍 Built-in search functionality
- 🌐 Custom domain support
- 📱 Mobile-responsive
- ⚡ CDN-powered, fast global delivery

**Setup Steps:**

1. **Create Account:**
   - Go to https://readthedocs.org
   - Sign in with your GitHub account

2. **Import Project:**
   - Click "Import a Project"
   - Select your `hpc-scaletest` repository
   - Click "Next"

3. **Configure:**
   - Read the Docs will automatically detect the `.readthedocs.yaml` file included in the package
   - It will build your documentation automatically

4. **Your Documentation URL:**
   ```
   https://hpc-scaletest.readthedocs.io
   ```

**That's it!** Read the Docs will rebuild your documentation automatically whenever you push changes.

---

### Option 2: GitHub Pages

**Best for:** Simple hosting, full control, static site

**Why Choose This:**
- 🆓 Free hosting
- 🔗 Easy to set up with GitHub Actions
- 🌐 Custom domain support
- 📦 Integrated with GitHub

**Setup Steps:**

**Automated with GitHub Actions (Recommended):**

This repository already includes a GitHub Actions workflow file at `.github/workflows/docs.yml` that automatically builds and deploys your documentation to GitHub Pages whenever changes are pushed to the `main` branch.

To enable GitHub Pages deployment:

1. Push this code to your GitHub repository
2. Go to your repository Settings → Pages
3. Under "Source", select "GitHub Actions"
4. The documentation will automatically be deployed on the next push to main

**Your Documentation URL:**
```
https://YOUR_USERNAME.github.io/hpc-scaletest
```

**Manual Deployment:**

If you prefer to deploy manually:

```bash
# 1. Extract the documentation package
unzip hpc-scaletest-documentation.zip
cd hpc-scaletest-documentation

# 2. Build the documentation
cd docs
pip install -r requirements.txt
make html

# 3. Deploy to GitHub Pages
cd build/html
git init
git add .
git commit -m "Deploy documentation"
git remote add origin https://github.com/YOUR_USERNAME/hpc-scaletest.git
git push -f origin HEAD:gh-pages

# 4. Enable GitHub Pages in repository settings
# Go to: Settings → Pages → Source: gh-pages branch
```
---

### Option 3: GitLab Pages

**Best for:** GitLab users, CI/CD integration

**Setup Steps:**

Create `.gitlab-ci.yml` in your repository root:

```yaml
pages:
  image: python:3.11
  stage: deploy
  script:
    - pip install -r docs/requirements.txt
    - cd docs && make html
    - mv build/html ../public
  artifacts:
    paths:
      - public
  only:
    - main
```

**Your Documentation URL:**
```
https://YOUR_USERNAME.gitlab.io/hpc-scaletest
```

---

### Option 4: Netlify

**Best for:** Modern hosting with deploy previews, great UX

**Why Choose This:**
- 🆓 Free tier available
- 👀 Deploy previews for pull requests
- ⚡ Fast CDN
- 🔐 HTTPS by default
- 🌐 Custom domains

**Setup Steps:**

1. Sign up at https://netlify.com
2. Click "New site from Git"
3. Connect your repository
4. Build settings:
   - **Build command:** `cd docs && make html`
   - **Publish directory:** `docs/build/html`
5. Click "Deploy site"

**Your Documentation URL:**
```
https://hpc-scaletest.netlify.app
```

(You can customize this to your own domain)

---

### Option 5: Self-Hosted

**Best for:** Internal documentation, enterprise environments

**Setup Steps:**

```bash
# 1. Build documentation
cd docs
pip install -r requirements.txt
make html

# 2. Copy to web server
scp -r build/html/* user@your-server:/var/www/html/hpc-scaletest-docs/

# Or use rsync
rsync -avz build/html/ user@your-server:/var/www/html/hpc-scaletest-docs/
```

**Nginx Configuration:**

```nginx
server {
    listen 80;
    server_name docs.your-domain.com;
    root /var/www/html/hpc-scaletest-docs;
    index index.html;
    
    location / {
        try_files $uri $uri/ =404;
    }
}
```

**Apache Configuration:**

```apache
<VirtualHost *:80>
    ServerName docs.your-domain.com
    DocumentRoot /var/www/html/hpc-scaletest-docs
    
    <Directory /var/www/html/hpc-scaletest-docs>
        Options Indexes FollowSymLinks
        AllowOverride All
        Require all granted
    </Directory>
</VirtualHost>
```

---

## 📝 Quick Start (Local Testing)

Before deploying, you can test locally:

```bash
# 1. Extract the package
unzip hpc-scaletest-documentation.zip
cd hpc-scaletest-documentation

# 2. Install dependencies
pip install -r docs/requirements.txt

# 3. Build documentation
cd docs
make html

# 4. Open in browser
# Open docs/build/html/index.html in your web browser
```

**Live Preview (Optional):**

```bash
cd docs
make livehtml
# Opens at http://127.0.0.1:8000
```

---

## 🎨 Customization Tips

### Change Theme

Edit `docs/source/conf.py`:

```python
html_theme = 'sphinx_rtd_theme'  # Current

# Other options:
# html_theme = 'alabaster'
# html_theme = 'sphinx_book_theme'
# html_theme = 'pydata_sphinx_theme'
```

### Add Your Logo

1. Place logo in `docs/source/_static/logo.png`
2. Edit `docs/source/conf.py`:

```python
html_logo = '_static/logo.png'
html_favicon = '_static/favicon.ico'
```

### Custom CSS

1. Create `docs/source/_static/custom.css`
2. Edit `docs/source/conf.py`:

```python
html_css_files = ['custom.css']
```

---

## 📂 What to Do Next

1. **Extract the documentation package** to your HPC-ScaleTest repository:
   ```bash
   cd /path/to/hpc-scaletest
   unzip hpc-scaletest-documentation.zip
   mv docs_package/docs .
   mv docs_package/.readthedocs.yaml .
   mv docs_package/README.md docs/
   ```

2. **Test locally:**
   ```bash
   cd docs
   pip install -r requirements.txt
   make html
   # Open build/html/index.html
   ```

3. **Choose hosting option** and deploy (I recommend Read the Docs)

4. **Update links** in your main README to point to your documentation

---

## 🆘 Troubleshooting

### Build Errors

```bash
# Clean and rebuild
cd docs
make clean
make html
```

### Missing Dependencies

```bash
pip install -r docs/requirements.txt
```

### Theme Not Found

```bash
pip install sphinx-rtd-theme
```

---

## 📚 Documentation Maintenance

- **Update regularly:** Keep docs in sync with code
- **Version tags:** Use semantic versioning (v1.0.0, v1.1.0)
- **Test builds:** Always build locally before pushing
- **Get feedback:** Ask users what needs clarification

---

## 🌟 Recommended: Read the Docs

**I strongly recommend Read the Docs** because:

1. Zero maintenance - automatic builds
2. Version management out of the box
3. Professional appearance
4. Trusted by major projects (Django, Python, Requests, etc.)
5. Free forever for open source
6. Great search functionality
7. Download options (PDF, ePub)

**It's literally 3 clicks to get started:**
1. Login to readthedocs.org with GitHub
2. Import your repository
3. Done!

---

## 📞 Support

If you need help:
- Check `DOCUMENTATION_GUIDE.md` in the package
- Read the Docs documentation: https://docs.readthedocs.io
- Sphinx documentation: https://www.sphinx-doc.org

---

**Your documentation is ready to go! Choose your hosting platform and deploy. Good luck! 🚀**
