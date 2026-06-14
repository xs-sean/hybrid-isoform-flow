# Public Release Checklist

Use a fresh clean repository or an orphan branch for the public archive. The current working tree excludes unpublished manuscript drafts, figures, tables, private outputs, and scratch scripts, but old Git commits may still contain those files.

Recommended release path:

```bash
git clone --no-local /path/to/idp_beta idp_beta_public
cd idp_beta_public
git checkout --orphan public-main
git add README.md LICENSE pyproject.toml .gitignore config container docs idp_beta results_example tests third_party workflow
git commit -m "Public IDP_BETA workflow release"
git remote add origin <new-public-repository-url>
git push -u origin public-main:main
```

Before pushing:

```bash
python3 -m unittest discover -s tests
python3 -m idp_beta.cli validate-config --no-check-files
python3 -m idp_beta.cli dry-run
find . -path ./.git -prune -o -type f \( -iname "*.docx" -o -iname "*.xlsx" -o -iname "*.pdf" -o -iname "*.tif" -o -iname "*.tiff" -o -iname "*.png" \) -print
rg -n "PRIVATE|UNPUBLISHED|LOCAL_ONLY|absolute/private/path/pattern" .
```

The `find` and `rg` commands should return no private manuscript/result files. The term `manuscript` may appear only in privacy/release documentation.
