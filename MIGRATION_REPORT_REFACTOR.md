# Migration Guide: sequana_report Extraction (Phase 2)

## Summary

HTML report infrastructure moved to `sequana_report` package. **No breaking changes** — old imports still work via backward compatibility layer.

## What Changed

### Dependency
Sequana now depends on `sequana_report>=0.1.0`. Installed automatically with `pip install sequana`.

### Import Paths

| Old (Deprecated) | New (Preferred) | Status |
|-----------------|-----------------|--------|
| `from sequana.modules_report import SequanaBaseModule` | `from sequana_report import SequanaBaseModule` | Works, shows DeprecationWarning |
| `from sequana.modules_report.base_module import ...` | `from sequana_report import SequanaBaseModule` | Works, shows DeprecationWarning |

### Report Modules
Modules like `rnadiff.py`, `coverage.py` stay in `sequana.modules_report` for now.

```python
# Still works (not moving yet):
from sequana.modules_report.summary import SequanaReport
from sequana.modules_report.coverage import CoverageModule
```

## No Action Required For...

- **Users:** No changes needed. Everything works as before.
- **Pipelines:** No changes needed. Backward compat handles it.
- **Tests:** Existing tests pass unchanged.

## Optional: Migrate Your Code

To avoid deprecation warnings, update imports:

### Before
```python
from sequana.modules_report import SequanaBaseModule
```

### After
```python
from sequana_report import SequanaBaseModule
```

That's it. Rest of your code stays the same.

## FAQ

**Q: Will old imports stop working?**
A: No. Compat layer stays for 1-2 releases minimum. Deprecation warnings warn of future removal.

**Q: Do I need to update my pipeline?**
A: No. Pipelines don't directly import from modules_report (they call CLI tools which handle it internally).

**Q: What if I use SequanaReport (not SequanaBaseModule)?**
A: SequanaReport stays in `sequana.modules_report.summary`. Keep importing from there.

**Q: How do I run report regeneration?**
A: Not available yet. Phase 4 feature. Check back in future releases.

## For Maintainers

### Running Tests
```bash
pytest test/modules_report/  # Uses compat layer, should all pass
```

### Checking Deprecation Warnings
```bash
python -W error::DeprecationWarning your_script.py  # Will fail on deprecated imports
```

### When to Deprecate Old Imports
- 1-2 releases after sequana_report 0.1.0 is stable
- Document removal in release notes
- Provide clear migration path

## Files Changed

- `requirements.txt` — added sequana_report>=0.1.0
- `sequana/modules_report/base_module.py` — now compat layer
- `sequana/modules_report/__init__.py` — re-exports from sequana_report

## See Also

- sequana_report REFACTOR.md — detailed architecture
- sequana_report GitHub: https://github.com/sequana/report
