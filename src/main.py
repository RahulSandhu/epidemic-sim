import sys

import matplotlib
import numpy as np
from PyQt6.QtWidgets import QApplication
from scipy.integrate import solve_ivp

matplotlib.use("QtAgg")

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QRadioButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)


# Define the SIR model equations
def sir_ode(t, y, beta, mu):
    S, I, R = y
    dSdt = -beta * S * I
    dIdt = beta * S * I - mu * I
    dRdt = mu * I
    return [dSdt, dIdt, dRdt]


# Define the SIS model equations
def sis_ode(t, y, beta, mu):
    S, I = y
    dSdt = mu * I - beta * S * I
    dIdt = beta * S * I - mu * I
    return [dSdt, dIdt]


# Define the SEIR model equations
def seir_ode(t, y, beta, alpha, mu):
    S, E, I, R = y
    dSdt = -beta * S * I
    dEdt = beta * S * I - alpha * E
    dIdt = alpha * E - mu * I
    dRdt = mu * I
    return [dSdt, dEdt, dIdt, dRdt]


# Define a custom Matplotlib Canvas widget
class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        fig.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.9)
        self.axes = fig.add_subplot(111)
        super().__init__(fig)


# Create the main application window class
class EpidemicApp(QMainWindow):
    def __init__(self):
        super().__init__()

        # Initialize values
        self.beta = 0.5
        self.mu = 0.1
        self.alpha = 0.2
        self.current_model = "SIR"
        self.scale_factor = 1000

        # Set up the main window
        self.setWindowTitle("Interactive Epidemic Modeler (SIR/SIS/SEIR)")
        self.resize(800, 750)

        # Create the main layout
        main_layout = QVBoxLayout()
        main_layout.setSpacing(0)
        main_layout.setContentsMargins(0, 0, 0, 0)

        # Create the Matplotlib canvas
        self.canvas = MplCanvas(self, width=5, height=4, dpi=100)
        main_layout.addWidget(self.canvas, 1)

        self.plot_lines = []
        self._create_tooltip()

        # Create the controls widget
        controls_widget = self._create_controls_widget()
        main_layout.addWidget(controls_widget, 0)

        # Set the central widget
        main_widget = QWidget()
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)

        # Connect signals and slots
        self._connect_signals()

        # Draw the initial plot
        self._update_plot()

    def _create_controls_widget(self):
        controls_container = QWidget()

        # Set white background for controls_container
        controls_container.setStyleSheet("background-color: white;")

        # This VBox will stack the grid on top of the radio buttons
        main_controls_vbox = QVBoxLayout()
        main_controls_vbox.setSpacing(10)
        main_controls_vbox.setContentsMargins(10, 5, 10, 10)

        # Grid layout for the sliders
        grid_layout = QGridLayout()
        grid_layout.setSpacing(10)

        large_font_style = "font-size: 14px;"
        slider_width = 300
        text_width = 70

        # Beta (Infection Rate)
        self.beta_label = QLabel("Beta (\u03b2):")
        self.beta_label.setStyleSheet(large_font_style)
        self.beta_slider = QSlider(Qt.Orientation.Horizontal)
        self.beta_slider.setRange(0, self.scale_factor)
        self.beta_slider.setValue(int(self.beta * self.scale_factor))
        self.beta_slider.setMinimumHeight(25)
        self.beta_slider.setFixedWidth(slider_width)
        self.beta_text = QLineEdit(f"{self.beta:.3f}")
        self.beta_text.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.beta_text.setFixedWidth(text_width)
        self.beta_text.setStyleSheet(large_font_style)

        grid_layout.addWidget(self.beta_label, 0, 0)
        grid_layout.addWidget(self.beta_slider, 0, 1)
        grid_layout.addWidget(self.beta_text, 0, 2)

        # Mu (Recovery Rate)
        self.mu_label = QLabel("Mu (\u03bc):")
        self.mu_label.setStyleSheet(large_font_style)
        self.mu_slider = QSlider(Qt.Orientation.Horizontal)
        self.mu_slider.setRange(0, self.scale_factor)
        self.mu_slider.setValue(int(self.mu * self.scale_factor))
        self.mu_slider.setMinimumHeight(25)
        self.mu_slider.setFixedWidth(slider_width)
        self.mu_text = QLineEdit(f"{self.mu:.3f}")
        self.mu_text.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.mu_text.setFixedWidth(text_width)
        self.mu_text.setStyleSheet(large_font_style)

        grid_layout.addWidget(self.mu_label, 1, 0)
        grid_layout.addWidget(self.mu_slider, 1, 1)
        grid_layout.addWidget(self.mu_text, 1, 2)

        # Alpha (Progression Rate)
        self.alpha_label = QLabel("Alpha (\u03b1):")
        self.alpha_label.setStyleSheet(large_font_style)
        self.alpha_slider = QSlider(Qt.Orientation.Horizontal)
        self.alpha_slider.setRange(0, self.scale_factor)
        self.alpha_slider.setValue(int(self.alpha * self.scale_factor))
        self.alpha_slider.setMinimumHeight(25)
        self.alpha_slider.setFixedWidth(slider_width)
        self.alpha_text = QLineEdit(f"{self.alpha:.3f}")
        self.alpha_text.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.alpha_text.setFixedWidth(text_width)
        self.alpha_text.setStyleSheet(large_font_style)

        # Position for the sliders
        grid_layout.addWidget(self.alpha_label, 2, 0)
        grid_layout.addWidget(self.alpha_slider, 2, 1)
        grid_layout.addWidget(self.alpha_text, 2, 2)

        # Add the grid of sliders to the vertical layout
        main_controls_vbox.addLayout(grid_layout)

        # Model Selector (Radio Buttons)
        self.sir_radio = QRadioButton("SIR")
        self.sir_radio.setStyleSheet(large_font_style)

        # Default selection
        self.sir_radio.setChecked(True)

        self.sis_radio = QRadioButton("SIS")
        self.sis_radio.setStyleSheet(large_font_style)

        self.seir_radio = QRadioButton("SEIR")
        self.seir_radio.setStyleSheet(large_font_style)

        # Create a HBox for the radio buttons
        model_hbox = QHBoxLayout()
        model_hbox.addStretch(2)
        model_hbox.addWidget(self.sir_radio)
        model_hbox.addStretch(1)
        model_hbox.addWidget(self.sis_radio)
        model_hbox.addStretch(1)
        model_hbox.addWidget(self.seir_radio)
        model_hbox.addStretch(2)

        # Add the radio button HBox to the VBox
        main_controls_vbox.addLayout(model_hbox)

        # This layout centers the entire VBox
        centering_layout = QHBoxLayout()
        centering_layout.addStretch(1)
        centering_layout.addLayout(main_controls_vbox)
        centering_layout.addStretch(1)
        controls_container.setLayout(centering_layout)

        return controls_container

    def _create_tooltip(self):
        """Creates a floating QLabel to act as a tooltip."""
        self.tooltip_label = QLabel(self.canvas)
        self.tooltip_label.setWindowFlags(Qt.WindowType.ToolTip)
        self.tooltip_label.setStyleSheet(
            f"""
            QLabel {{
                background-color: white;
                color: black;
                border: 1px solid black;
                padding: 4px;
                border-radius: 3px;
                font-size: 13px;
            }}
        """
        )
        self.tooltip_label.hide()
        self.tooltip_label.setAlignment(Qt.AlignmentFlag.AlignLeft)

    def _on_hover(self, event):
        """Handles mouse motion over the Matplotlib canvas to show a tooltip."""
        if not event.inaxes or event.guiEvent is None:
            self.tooltip_label.hide()
            return

        if not self.plot_lines:
            self.tooltip_label.hide()
            return

        min_dist = float("inf")
        closest_point = None
        line_label = ""

        for line in self.plot_lines:
            xdata, ydata = line.get_data()
            if len(xdata) == 0:
                continue

            idx = np.argmin(np.abs(xdata - event.xdata))
            x, y = xdata[idx], ydata[idx]

            display_x, display_y = self.canvas.axes.transData.transform((x, y))
            dist = np.sqrt((display_x - event.x) ** 2 + (display_y - event.y) ** 2)

            if dist < min_dist:
                min_dist = dist
                closest_point = (x, y)
                line_label = line.get_label()

        # Only show the tooltip if the mouse is close to a data point
        if min_dist < 20:
            x, y = closest_point
            text = f"{line_label}\nDay: {x:.1f}\nValue: {y:.3f}"
            self.tooltip_label.setText(text)

            pos = event.guiEvent.globalPosition().toPoint()

            self.tooltip_label.move(pos.x() + 10, pos.y() + 10)
            self.tooltip_label.show()
            self.tooltip_label.adjustSize()
        else:
            self.tooltip_label.hide()

    def _connect_signals(self):
        """Connect all signals to their handler functions (slots)."""
        self.beta_slider.valueChanged.connect(self._beta_slider_changed)
        self.beta_text.editingFinished.connect(self._beta_text_changed)

        self.mu_slider.valueChanged.connect(self._mu_slider_changed)
        self.mu_text.editingFinished.connect(self._mu_text_changed)

        self.alpha_slider.valueChanged.connect(self._alpha_slider_changed)
        self.alpha_text.editingFinished.connect(self._alpha_text_changed)

        self.canvas.mpl_connect("motion_notify_event", self._on_hover)

        # Connect the new radio buttons
        self.sir_radio.toggled.connect(self._on_model_selected)
        self.sis_radio.toggled.connect(self._on_model_selected)
        self.seir_radio.toggled.connect(self._on_model_selected)

    # New slot for Radio Buttons
    def _on_model_selected(self):
        """Handles when a radio button is toggled."""
        # Check which button is currently checked
        if self.sir_radio.isChecked():
            new_model = "SIR"
        elif self.sis_radio.isChecked():
            new_model = "SIS"
        elif self.seir_radio.isChecked():
            new_model = "SEIR"
        else:
            return

        # Only update if the model actually changed
        if self.current_model != new_model:
            self.current_model = new_model
            self._update_plot()

    # Slot functions for Beta

    def _beta_slider_changed(self, value):
        self.beta = value / self.scale_factor
        self.beta_text.blockSignals(True)
        self.beta_text.setText(f"{self.beta:.3f}")
        self.beta_text.blockSignals(False)
        self._update_plot()

    def _beta_text_changed(self):
        try:
            val = float(self.beta_text.text())
        except ValueError:
            val = self.beta
        val = max(0.0, min(1.0, val))
        self.beta = val
        self.beta_text.setText(f"{self.beta:.3f}")
        self.beta_slider.blockSignals(True)
        self.beta_slider.setValue(int(self.beta * self.scale_factor))
        self.beta_slider.blockSignals(False)
        self._update_plot()

    # Slot functions for Mu

    def _mu_slider_changed(self, value):
        self.mu = value / self.scale_factor
        self.mu_text.blockSignals(True)
        self.mu_text.setText(f"{self.mu:.3f}")
        self.mu_text.blockSignals(False)
        self._update_plot()

    def _mu_text_changed(self):
        try:
            val = float(self.mu_text.text())
        except ValueError:
            val = self.mu
        val = max(0.0, min(1.0, val))
        self.mu = val
        self.mu_text.setText(f"{self.mu:.3f}")
        self.mu_slider.blockSignals(True)
        self.mu_slider.setValue(int(self.mu * self.scale_factor))
        self.mu_slider.blockSignals(False)
        self._update_plot()

    # Slot functions for Alpha

    def _alpha_slider_changed(self, value):
        self.alpha = value / self.scale_factor
        self.alpha_text.blockSignals(True)
        self.alpha_text.setText(f"{self.alpha:.3f}")
        self.alpha_text.blockSignals(False)
        self._update_plot()

    def _alpha_text_changed(self):
        try:
            val = float(self.alpha_text.text())
        except ValueError:
            val = self.alpha
        val = max(0.0, min(1.0, val))
        self.alpha = val
        self.alpha_text.setText(f"{self.alpha:.3f}")
        self.alpha_slider.blockSignals(True)
        self.alpha_slider.setValue(int(self.alpha * self.scale_factor))
        self.alpha_slider.blockSignals(False)
        self._update_plot()

    # Plotting Function

    def _update_plot(self):
        t_span = [0, 160]
        t_points = np.linspace(t_span[0], t_span[1], 300)

        I0 = 0.01
        R0 = 0.0
        E0 = 0.01

        if self.current_model == "SIR":
            S0 = 1.0 - I0 - R0
            y0 = [S0, I0, R0]
            ode_func = sir_ode
            args = (self.beta, self.mu)

        elif self.current_model == "SIS":
            S0 = 1.0 - I0
            y0 = [S0, I0]
            ode_func = sis_ode
            args = (self.beta, self.mu)

        elif self.current_model == "SEIR":
            S0 = 1.0 - I0 - R0 - E0
            y0 = [S0, E0, I0, R0]
            ode_func = seir_ode
            args = (self.beta, self.alpha, self.mu)

        sol = solve_ivp(ode_func, t_span, y0, t_eval=t_points, args=args)

        t = sol.t
        S, I, R, E = (np.zeros_like(t) for _ in range(4))

        if self.current_model == "SIR":
            S, I, R = sol.y
        elif self.current_model == "SIS":
            S, I = sol.y
        elif self.current_model == "SEIR":
            S, E, I, R = sol.y

        self.canvas.axes.clear()

        self.plot_lines.clear()

        (line_s,) = self.canvas.axes.plot(t, S, label="Susceptible (S)", color="green")
        self.plot_lines.append(line_s)

        (line_i,) = self.canvas.axes.plot(t, I, label="Infectious (I)", color="red")
        self.plot_lines.append(line_i)

        if self.current_model == "SEIR":
            (line_e,) = self.canvas.axes.plot(t, E, label="Exposed (E)", color="orange")
            self.plot_lines.append(line_e)

        if self.current_model in ["SIR", "SEIR"]:
            (line_r,) = self.canvas.axes.plot(t, R, label="Recovered (R)", color="blue")
            self.plot_lines.append(line_r)

        R0_val = self.beta / self.mu if self.mu > 0 else float("inf")

        title = f"Interactive {self.current_model} Model (R\u2080 \u2248 {R0_val:.2f})"
        self.canvas.axes.set_title(title)
        self.canvas.axes.set_xlabel("Time (Days)")
        self.canvas.axes.set_ylabel("Fraction of Population")
        self.canvas.axes.legend(loc="upper right")
        self.canvas.axes.grid(True)
        self.canvas.axes.set_ylim(-0.01, 1.01)

        self.canvas.draw()


def main():
    """Main function to run the application."""
    app = QApplication(sys.argv)
    window = EpidemicApp()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
