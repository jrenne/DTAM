# Security policy

Please do not place passwords, API keys, tokens, or licensed data in DTAM
source files, examples, issues, or pull requests.

Data-construction scripts must obtain credentials from environment variables.
For FRED, use `FRED_API_KEY`; each user should supply their own key.

If a credential is accidentally committed, remove it from the current source
and revoke or rotate it immediately. Removing a value from the latest commit
does not remove it from Git history.

Report a suspected credential exposure privately to the package maintainer at
<jean-paul.renne@unil.ch>. Do not include the credential itself in the report.
