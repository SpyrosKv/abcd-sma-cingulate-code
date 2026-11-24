"""Check the code-only release tree (or staged Git files) without printing matched values."""
from pathlib import Path
import argparse, re, subprocess, sys

def audit(root, paths):
    allowed = {'.R', '.md', '.Rproj', '.cff', '.py', '.lock'}
    permitted_csv = {'docs/source_provenance.csv', 'docs/package_manifest.csv', 'docs/validated_package_versions.csv'}
    findings = []
    patterns = {
        'participant identifier': r'(?:sub-|NDAR_INV)[A-Z0-9]{6,}',
        'GitHub credential': r'gh[pousr]_[A-Za-z0-9]{20,}',
        'private key': r'-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----',
        'view-only access token': r'view_only=[a-fA-F0-9]{20,}',
        'personal Windows path': r'[A-Za-z]:[\\/](?:Users|Documents and Settings)[\\/]',
        'cloud access key': r'AKIA[A-Z0-9]{16}',
    }
    for path in paths:
        rel = path.relative_to(root).as_posix()
        if path.is_symlink():
            findings.append((rel, 'symlink is not allowed')); continue
        if path.name in {'.gitignore', '.gitattributes'}: continue
        if ((path.name.startswith('config') and path.suffix == '.R' and path.name != 'config.example.R') or path.name in {'.Renviron', '.Rhistory', '.RData'}):
            findings.append((rel, 'local/private file')); continue
        if path.suffix not in allowed and rel not in permitted_csv and rel != 'LICENSE':
            findings.append((rel, 'file type not in code-only allowlist')); continue
        try: text = path.read_text(encoding='utf-8')
        except (UnicodeError, OSError):
            findings.append((rel, 'not readable UTF-8 source text')); continue
        for label, pattern in patterns.items():
            if re.search(pattern, text): findings.append((rel, label))
    return findings

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--staged', action='store_true')
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    if args.staged:
        result = subprocess.run(['git', 'ls-files', '--cached', '-z'], cwd=root,
                                check=True, capture_output=True)
        paths = [root / x.decode() for x in result.stdout.split(b'\0') if x]
        # Also reject staged bytes differing from the inspected working copy.
        dirty = subprocess.run(['git', 'diff', '--name-only', '-z'], cwd=root,
                               check=True, capture_output=True)
        if dirty.stdout:
            print('FAIL: staged/working files differ; inspect the staged content directly.'); return 1
    else:
        paths = [p for p in root.rglob('*') if p.is_file() and '.git' not in p.parts
                 and '.Rproj.user' not in p.parts and '__pycache__' not in p.parts]
    findings = audit(root, paths)
    for path, reason in findings: print(f'FAIL: {path}: {reason}')
    if findings: return 1
    print(f'PASS: {len(paths)} code/documentation files inspected. Manual review still required.')
    return 0

if __name__ == '__main__': sys.exit(main())
