#!/usr/bin/env python3
import json
import sys

def print_tree(data, indent=0):
    if isinstance(data, dict):
        for key, value in data.items():
            print("  " * indent + f"{key}:")
            print_tree(value, indent + 1)
    elif isinstance(data, list):
        for i, value in enumerate(data):
            print("  " * indent + f"[{i}]:")
            print_tree(value, indent + 1)
    else:
        print("  " * indent + str(data))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <json_file>")
        sys.exit(1)
    
    json_file = sys.argv[1]
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    print_tree(data)
