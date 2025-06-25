attr_dict = {
    "u": {
        "standard_name": "eastward_wind",
        "long_name": "Eastward wind component",
        "units": "m s-1",
    },
    "v": {
        "standard_name": "northward_wind",
        "long_name": "Northward wind component",
        "units": "m s-1",
    },
    "q": {
        "standard_name": "specific_humidity",
        "long_name": "Specific humidity",
        "units": "kg kg-1",
    },
    "ta": {
        "standard_name": "air_temperature",
        "long_name": "Air temperature",
        "units": "K",
    },
    "p": {
        "standard_name": "air_pressure",
        "long_name": "Atmospheric pressure",
        "units": "Pa",
    },
    "x": {
        "long_name": "x",
        "units": "m",
        "description": "Distance of sonde longitude to mean circle longitude",
    },
    "y": {
        "long_name": "y",
        "units": "m",
        "description": "Distance of sonde latitude to mean circle latitude",
    },
    "u_mean": {
        "long_name": "circle mean of u component of winds",
        "units": "m s-1",
    },
    "v_mean": {
        "long_name": "circle mean of v component of winds",
        "units": "m s-1",
    },
    "q_mean": {
        "long_name": "circle mean of specific humidity",
        "units": "kg kg-1",
    },
    "ta_mean": {
        "long_name": "circle mean of air temperature",
        "units": "K",
    },
    "vor": {
        "long_name": "Area-averaged horizontal relative vorticity",
        "units": "s-1",
        "standard_name": "atmosphere_relative_vorticity",
    },
    "div": {
        "long_name": "Area-averaged horizontal mass divergence",
        "units": "s-1",
        "standard_name": "divergence_of_wind",
    },
    "omega": {
        "long_name": "Area-averaged atmospheric pressure velocity (omega)",
        "units": "hPa hr-1",
        "standard_name": "vertical_air_velocity_expressed_as_tendency_of_pressure",
    },
    "analysis_type": {
        "long_name": "(IFS) Analysis type",
        "description": "orc: ORCESTRA measurements, obs: IFS observations, an: IFS analysis, fg: IFS first guess",
    },
    "circle_id": {
        "long_name": "circle identifier",
        "description": "Unique circle ID from flight segmentation",
    },
}

circle_attr_dict = {
    f"{var}_d{var}dx": {
        "long_name": f"zonal gradient of {attr_dict[var]['long_name']}",
        "units": f"{attr_dict[var]['units']} m-1",
        "standard_name": f"derivative_of_{attr_dict[var]['standard_name']}_wrt_x",
    }
    for var in ["u", "v", "q", "ta"]
}
circle_attr_dict.update(
    {
        f"{var}_d{var}dy": {
            "long_name": f"meridional gradient of {attr_dict[var]['long_name']}",
            "units": f"{attr_dict[var]['units']} m-1",
            "standard_name": f"derivative_of_{attr_dict[var]['standard_name']}_wrt_y",
        }
        for var in ["u", "v", "q", "ta"]
    }
)
circle_attr_dict.update(
    {
        f"{var}_mean": {
            "long_name": f"circle mean of {attr_dict[var]['long_name']}",
            "units": attr_dict[var]["units"],
        }
        for var in ["u", "v", "q", "ta"]
    }
)
