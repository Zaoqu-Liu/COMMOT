# Bug Fix: httr2 Compatibility

## Issue
Fixed compatibility issue with httr2 >= 1.2.0 that caused download failures with error:
```
unused arguments (req_prep, resend_count = n)
```

## Solution
Updated `R/httr2_keep_handle.R` to dynamically detect and adapt to httr2 API changes. The patch function now correctly handles both old (3 parameters) and new (5 parameters) versions of `req_perform1`.

## Changes
- Modified `R/httr2_keep_handle.R`: Dynamic function signature detection
- Enhanced `R/internals.R`: Improved error handling

## Testing
All functions now work correctly with httr2 1.2.2:
- Ensembl organism data download ✓
- DoRothEA transcription factor data ✓
- OmniPath interactions ✓
- Protein annotations ✓

## Maintainer
Zaoqu Liu (liuzaoqu@163.com)

