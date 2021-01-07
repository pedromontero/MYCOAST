# MYCOAST
INTECMAR code developed for MYCOAST project.

SCRIPTS:
-------

**nc_struct_to_json:** Read a netcdf file and create a json file, 
describing the structure of it.

**json_struct_to_nc:** Read a json file with the structure for  a
netcdf and create empty netcf file with this structure.

The aim of these complementary scripts is the next:

In order to write a netcdf file with the specifications of MyCOAST, the 
structure of this file (metadata and variables) is stored in a json file:
```json 
{
    "dimensions": 
    [
        {"name": "TIME","size": 24 },
        {"name": "FREQUENCY","size": 14 }
    ],
    "global attributes":
    [
        {"name": "platform_code","value": "6100430","type": "str"}
        {"name": "platform_name","value": "Dragonera buoy","type": "str"}
    ],
    "variables":
     [
        {"name": "TIME","type": "f8","fill_value": null,
         "dimensions": ["TIME"],
         "attributes": 
             [
                 {"name": "long_name", "value": "Time","type": "str"},
                 {"name": "standard_name","value": "time","type": "str"},
                 {"name": "units","value": "days since 1950-01-01T00:00:00Z","type"}
            ]
        }      
    ]
}
```
Then this structure is read to create other like-structure netCDF files and filling them
with actual data. Little changes in metadata can be changed by hand or by code.

In order to create the first like-structure netcdf file, nc_struct_to_json reads a 
current netcdf file and reproduces its structure to the json structure file.

   


