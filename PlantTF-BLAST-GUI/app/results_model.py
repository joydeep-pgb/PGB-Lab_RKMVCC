"""
results_model.py
------------------
Qt table model for the results table (one row per query sequence) and a
QSortFilterProxyModel subclass that supports free-text search, a TF-family
filter, and a "matched only" toggle.
"""

from typing import List, Optional

from PyQt6.QtCore import QAbstractTableModel, QModelIndex, QSortFilterProxyModel, Qt

from .blast_engine import QueryResult

COLUMNS = [
    ("Query ID", "query_id"),
    ("TF Family", "family"),
    ("Best Hit ID", "hit_id"),
    ("E-value", "evalue"),
    ("Identity (%)", "pident"),
    ("Query Cov (%)", "qcov"),
    ("Description", "description"),
]


def _fmt_evalue(v) -> str:
    if v is None:
        return "-"
    if v == 0:
        return "0.0"
    return f"{v:.2e}"


def _fmt_pct(v) -> str:
    if v is None:
        return "-"
    return f"{v:.1f}"


class ResultsTableModel(QAbstractTableModel):
    def __init__(self, results: Optional[List[QueryResult]] = None, parent=None):
        super().__init__(parent)
        self._results: List[QueryResult] = results or []

    # -- data loading -----------------------------------------------------
    def set_results(self, results: List[QueryResult]):
        self.beginResetModel()
        self._results = results
        self.endResetModel()

    def result_at(self, row: int) -> Optional[QueryResult]:
        if 0 <= row < len(self._results):
            return self._results[row]
        return None

    # -- required overrides ------------------------------------------------
    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self._results)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(COLUMNS)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            return COLUMNS[section][0]
        return str(section + 1)

    def data(self, index: QModelIndex, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        r = self._results[index.row()]
        col_key = COLUMNS[index.column()][1]

        if role == Qt.ItemDataRole.DisplayRole:
            if not r.matched:
                if col_key == "query_id":
                    return r.query_id
                if col_key == "family":
                    return "No hit"
                return "-"
            value = getattr(r, col_key)
            if col_key == "evalue":
                return _fmt_evalue(value)
            if col_key in ("pident", "qcov"):
                return _fmt_pct(value)
            return value if value not in (None, "") else "-"

        if role == Qt.ItemDataRole.TextAlignmentRole:
            if col_key in ("evalue", "pident", "qcov"):
                return Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
            return Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter

        if role == Qt.ItemDataRole.ForegroundRole and not r.matched:
            from PyQt6.QtGui import QColor

            return QColor("#999999")

        # raw sort role for numeric columns
        if role == Qt.ItemDataRole.UserRole:
            if not r.matched:
                return None
            value = getattr(r, col_key)
            return value

        return None


class ResultsFilterProxyModel(QSortFilterProxyModel):
    """Supports free-text search across all columns, a family filter and a
    'matched only' toggle."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._search_text = ""
        self._family_filter = "All families"
        self._matched_only = False
        self.setSortRole(Qt.ItemDataRole.DisplayRole)

    def set_search_text(self, text: str):
        self._search_text = text.strip().lower()
        self.invalidateFilter()

    def set_family_filter(self, family: str):
        self._family_filter = family
        self.invalidateFilter()

    def set_matched_only(self, matched_only: bool):
        self._matched_only = matched_only
        self.invalidateFilter()

    def filterAcceptsRow(self, source_row, source_parent):
        model: ResultsTableModel = self.sourceModel()
        r = model.result_at(source_row)
        if r is None:
            return False

        if self._matched_only and not r.matched:
            return False

        if self._family_filter not in ("All families", "", None):
            if r.family != self._family_filter:
                return False

        if self._search_text:
            haystack = " ".join(
                str(x)
                for x in (
                    r.query_id,
                    r.family,
                    r.hit_id,
                    r.description,
                    r.species,
                )
            ).lower()
            if self._search_text not in haystack:
                return False

        return True

    # numeric-aware sorting for E-value / identity / coverage columns
    def lessThan(self, left, right):
        col_key = COLUMNS[left.column()][1]
        if col_key in ("evalue", "pident", "qcov"):
            model: ResultsTableModel = self.sourceModel()
            lval = model.data(left, Qt.ItemDataRole.UserRole)
            rval = model.data(right, Qt.ItemDataRole.UserRole)
            if lval is None and rval is None:
                return False
            if lval is None:
                return False
            if rval is None:
                return True
            return lval < rval
        return super().lessThan(left, right)
