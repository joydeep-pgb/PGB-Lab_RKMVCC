"""
family_chart.py
-----------------
A small, dependency-free horizontal bar chart widget (drawn with QPainter)
that shows how many query sequences were classified into each TF family.
Avoids requiring matplotlib/QtCharts just for this one summary view.
"""

from typing import List, Tuple

from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtGui import QColor, QFont, QPainter
from PyQt6.QtWidgets import QScrollArea, QSizePolicy, QWidget

_PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
]


class FamilyBarChart(QWidget):
    ROW_HEIGHT = 22
    LABEL_WIDTH = 130
    COUNT_WIDTH = 40
    MARGIN = 8

    def __init__(self, parent=None):
        super().__init__(parent)
        self._data: List[Tuple[str, int]] = []
        self.setMinimumHeight(60)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

    def set_data(self, data: List[Tuple[str, int]]):
        self._data = data
        total_height = max(60, self.MARGIN * 2 + self.ROW_HEIGHT * max(1, len(data)))
        self.setMinimumHeight(total_height)
        self.resize(self.width(), total_height)
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        if not self._data:
            painter.setPen(QColor("#999999"))
            painter.drawText(
                self.rect(),
                Qt.AlignmentFlag.AlignCenter,
                "No TF family data yet. Run a BLAST search to see the distribution.",
            )
            return

        max_count = max(c for _, c in self._data) or 1
        bar_area_x0 = self.MARGIN + self.LABEL_WIDTH
        bar_area_w = max(20, self.width() - bar_area_x0 - self.COUNT_WIDTH - self.MARGIN)

        font = QFont()
        font.setPointSize(9)
        painter.setFont(font)

        for i, (family, count) in enumerate(self._data):
            y = self.MARGIN + i * self.ROW_HEIGHT
            color = QColor(_PALETTE[i % len(_PALETTE)])

            # label
            painter.setPen(QColor("#333333"))
            label_rect = QRectF(self.MARGIN, y, self.LABEL_WIDTH - 6, self.ROW_HEIGHT - 4)
            metrics = painter.fontMetrics()
            elided = metrics.elidedText(
                family, Qt.TextElideMode.ElideRight, int(label_rect.width())
            )
            painter.drawText(label_rect, Qt.AlignmentFlag.AlignVCenter, elided)

            # bar
            bar_w = max(2, int(bar_area_w * count / max_count))
            bar_rect = QRectF(bar_area_x0, y + 3, bar_w, self.ROW_HEIGHT - 8)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(color)
            painter.drawRoundedRect(bar_rect, 3, 3)

            # count
            painter.setPen(QColor("#333333"))
            count_rect = QRectF(
                bar_area_x0 + bar_w + 6, y, self.COUNT_WIDTH, self.ROW_HEIGHT - 4
            )
            painter.drawText(
                count_rect, Qt.AlignmentFlag.AlignVCenter, str(count)
            )


class FamilyChartScrollArea(QScrollArea):
    """Scrollable container for FamilyBarChart, since there can be 50+ families."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.chart = FamilyBarChart()
        self.setWidget(self.chart)
        self.setWidgetResizable(True)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

    def set_data(self, data: List[Tuple[str, int]]):
        self.chart.set_data(data)
