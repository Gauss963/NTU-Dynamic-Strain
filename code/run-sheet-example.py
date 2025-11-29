import json
import datetime
import numpy as np

serial_number = "t9999"
temperature = 25
humidity = 60

runsheet = {
    "date": str(datetime.date.today()),
    "serial_number": serial_number,
    "experimenter": ["Alice", "Bob", "Charlie", "David"],
    "temperature": temperature,
    "humidity": humidity,
    "sample_id": "PMMA-50mm",
    "channels": {
        "ch_a1": {
            "sensor-type": "PZT",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "z",
            "location": "stationary-block",
            "coordinates": [22.0, 10.0, 0.0],
        },
        "ch_a2": {
            "sensor-type": "PZT",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "z",
            "location": "stationary-block",
            "coordinates": [22.0, 20.0, 0.0],
        },
        "ch_a3": {
            "sensor-type": "PZT",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "z",
            "location": "stationary-block",
            "coordinates": [22.0, 30.0, 0.0],
        },
        "ch_a4": {
            "sensor-type": "PZT",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "z",
            "location": "stationary-block",
            "coordinates": [22.0, 40.0, 0.0],
        },
        "ch_a5": {
            "sensor-type": "ACC",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "x",
            "location": "moving-block",
            "coordinates": [20.0, 30.0, 0.0],
        },
        "ch_a6": {
            "sensor-type": "ACC",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "y",
            "location": "moving-block",
            "coordinates": [20.0, 30.0, 0.0],
        },
        "ch_a7": {
            "sensor-type": "ACC",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "z",
            "location": "moving-block",
            "coordinates": [20.0, 30.0, 0.0],
        },
        "ch_b1": {
            "sensor-type": "pressure",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "normal",
            "location": "moving-block",
            "coordinates": [0.0, 0.0, 0.0],
        },
        "ch_b2": {
            "sensor-type": "pressure",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "shear",
            "location": "moving-block",
            "coordinates": [0.0, 0.0, 0.0],
        },
        "ch_b3": {
            "sensor-type": "LVDT",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "y",
            "location": "moving-block",
            "coordinates": [10.0, 0.0, 0.0],
        },
        "ch_b4": {
            "sensor-type": "eddy",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "y",
            "location": "moving-block",
            "coordinates": [20.0, 5.0, 0.0],
        },
        "ch_b5": {
            "sensor-type": "eddy",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "y",
            "location": "moving-block",
            "coordinates": [20.0, 15.0, 0.0],
        },
        "ch_b6": {
            "sensor-type": "eddy",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "y",
            "location": "moving-block",
            "coordinates": [20.0, 25.0, 0.0],
        },
        "ch_b7": {
            "sensor-type": "eddy",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "y",
            "location": "moving-block",
            "coordinates": [20.0, 35.0, 0.0],
        },
        "ch_b8": {
            "sensor-type": "eddy",
            "range": "500mV",
            "hysteresis": "0.001V",
            "direction": "y",
            "location": "moving-block",
            "coordinates": [20.0, 45.0, 0.0],
        }
    },
    "notes": "This is an example experiment run sheet.",
    "cycle": ["32MPa", "15MPa", "8MPa", "4MPa"]
}