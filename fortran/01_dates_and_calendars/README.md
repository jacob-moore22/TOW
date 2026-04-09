# Chapter 1: Dates and Calendars

Calendar computations and date conversion routines.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `flmoon` | 1.0 | Calculate phases of the moon by date |
| `julday` | 1.1 | Julian Day number from calendar date |
| `badluk` | 1.1 | Find Friday the 13th when the moon is full (standalone program) |
| `caldat` | 1.1 | Calendar date from Julian day number |

## Notes

- `badluk` is a standalone `PROGRAM` (not a subroutine) that uses `julday` and `flmoon` to search for full moons on Friday the 13th between 1900 and 2000.
- `julday` and `caldat` are inverse operations: Julian day to/from calendar date.
- `flmoon` computes the Julian date and fractional time of a given lunar phase.
