## Def convert help
usage_mongo="""
gff-toolbox:

            Ask mongo

This command allows you to ask/interrogate the mongo db created from your GFF annotation file

usage:
    gff-toolbox ask-mongo [ -h|--help ] [ --list_dbs ] [ --list_collections <db_name> ]

options:

    -h, --help                                               Show this screen
    --list_dbs                                               List available mongo dbs in your system
    --list_collections=<db_name>                             Check the available collections in a given database

example:
"""

##################################
### Loading Necessary Packages ###
##################################
import sys
import os
import re
import urllib.request, urllib.parse, urllib.error
import json
import pymongo
from pymongo import MongoClient
import pathlib

########################
### Useful functions ###
########################

#########################
### Check your mongos ###
#########################
def check_mongos():

    # Start message
    print("""
    The available mondo dbs found in your system are:\n
    """)

    # Print databases
    os.system('mongo --quiet --eval  "printjson(db.adminCommand(\'listDatabases\'))"')

#######################################
### Check collections in a database ###
#######################################
def check_db(db_name):

    # Create connection
    client = MongoClient()

    # Open Database
    db = client[db_name]

    # Check available collections
    cols = db.list_collection_names()
    print(f"\nAll the available collections found in the {db_name} database are given in the list: {cols}\n")

################
### Def main ###
################
def ask_mongo(list_dbs, db_name):

    if list_dbs:
        check_mongos()

    elif db_name != None:
        check_db(db_name)
