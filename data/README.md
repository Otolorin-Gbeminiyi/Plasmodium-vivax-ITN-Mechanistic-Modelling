# Input data (not included)

The original trial screening data are deliberately **not packaged** with the public code. Public release requires permission from the relevant data owners and must meet applicable ethics, institutional and data-sharing conditions.

Two CSV files are read interactively by the single R script:

1. **CONTROL** arm CSV — chosen first.
2. **INTERVENTION** arm CSV — chosen second.

Each CSV requires these case-sensitive columns:

- `Date_of_the_survey` — a survey date; ISO `YYYY-MM-DD` is preferred.
- `number_tested` — numeric screening denominator.
- `number_positive` — numeric count of positives, between zero and `number_tested`.

An optional `cum_positive` column may be present. The script groups records by `Date_of_the_survey` unless configured to use a separate survey-round column. The source currently treats all October screening dates as passive case finding, excludes them from the primary likelihood and retains them descriptively. Primary-fitting figures require the year 2001.

**Data access:** The repository does not grant rights to access, use or redistribute the underlying trial data. Contact the appropriate research/data custodian.

**Safety note:** `.gitignore` prevents routine Git tracking of certain file formats; it cannot stop someone from deliberately uploading a file through the GitHub website. Always check the files before publishing.
