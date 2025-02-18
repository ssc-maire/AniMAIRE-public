import pytest
import pandas as pd
import numpy as np
from AniMAIRE.dose_plotting import plot_dose_map_contours, create_single_dose_map_plot_plt, plot_dose_map

@pytest.fixture
def sample_dose_map():
    data = {
        "longitude": np.linspace(0, 360, 10),
        "longitudeTranslated": np.linspace(-180, 180, 10),
        "latitude": np.linspace(-90, 90, 10),
        "altitude (km)": np.full(10, 0.0),
        "edose": np.random.rand(10),
        "SEU": np.random.rand(10),
        "SEL": np.random.rand(10),
    }
    return pd.DataFrame(data)

def test_plot_dose_map_contours(sample_dose_map):
    plot_dose_map_contours(sample_dose_map, levels=3)

def test_create_single_dose_map_plot_plt(sample_dose_map):
    create_single_dose_map_plot_plt(sample_dose_map)

def test_plot_dose_map(sample_dose_map):
    plot_dose_map(sample_dose_map, plot_title="Test Plot", plot_contours=True, levels=3)
