# Account Quick Reference

## 🎯 How to Specify Your Account

### ✅ Method 1: Command Line (RECOMMENDED)

Use the `--account` flag:

```bash
python hpc_auto.py REPO_URL \
    --scaling weak \
    --nodes 4 \
    --partition-name dcgp \
    --account YOUR_ACCOUNT_NAME \  # ← Your account here
    --verbose
```

### 📋 Common Examples

```bash
# CINECA staff
--account cin_staff

# Project account
--account IscrC_MyProject

# Budget account  
--account cin0123

# Personal account
--account username_budget
```

## 🔍 Check Your Accounts

```bash
# List available accounts
saldo -b

# Check account usage
saldo -r

# View past jobs with accounts
sacct --format=Account,JobID,JobName,State | head -20
```

## ✅ Your Current Command

You're already doing it correctly:

```bash
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git \
    --scaling weak \
    --nodes 16 \
    --partition-name dcgp \
    --account cin_staff \  # ✓ Your account specified here
    --modules "hdf5/1.14.3--intel-oneapi-mpi--2021.12.1--oneapi--2024.1.0" \
    --input-file inputfiles/os-stdin \
    --verbose
```

## 📝 Verification

After running, check the generated job script:

```bash
cat output/iPIC3D-CPU-NS_weak_*/nodes_1/job.sh | grep "#SBATCH -A"
```

Should show:
```bash
#SBATCH -A cin_staff  # ✓ Correct
```

NOT:
```bash
#SBATCH -A cin_X  # ✗ Default placeholder
```

## 🎓 Account Priority

1. **`--account` flag** ← Highest priority (use this!)
2. System config default (`cin_staff` for DCGP)
3. Built-in default (`cin_X` - placeholder only)

## 💡 Pro Tip

**Always use `--account YOUR_ACCOUNT`** to be explicit and avoid confusion!

## 🚀 Full Command Template

```bash
python hpc_auto.py <GITHUB_REPO_URL> \
    --scaling <weak|strong> \
    --nodes <MAX_NODES> \
    --partition-name <dcgp|booster> \
    --account <YOUR_ACCOUNT> \
    --modules "<MODULE1> <MODULE2>" \
    --input-file <PATH_TO_INPUT> \
    --verbose
```

Replace `<...>` with your values!

## 📞 Need Help?

Check: `python hpc_auto.py --help | grep -A 2 account`

Output:
```
--account ACCOUNT   Account/project name (default: cin_X)
```

**You can use any valid account name!** 🎉
