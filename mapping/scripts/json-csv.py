import json
import csv

def flatten_json(nested_json, parent_key='', sep='.'):
    """
    Flattens a nested JSON object.
    :param nested_json: The nested JSON object.
    :param parent_key: String representing the parent key.
    :param sep: Separator for the keys in the flattened structure.
    :return: A flattened dictionary.
    """
    items = {}
    for k, v in nested_json.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.update(flatten_json(v, new_key, sep=sep))
        elif isinstance(v, list):
            for i, item in enumerate(v):
                if isinstance(item, dict):
                    items.update(flatten_json(item, f"{new_key}{sep}{i}", sep=sep))
                else:
                    items[f"{new_key}{sep}{i}"] = item
        else:
            items[new_key] = v
    return items

def json_to_csv(json_file, csv_file):
    """
    Convert the JSON file to a CSV file.
    :param json_file: The input JSON file.
    :param csv_file: The output CSV file.
    """
    # Open and load the JSON data
    with open(json_file, "r") as f:
        # Since the file is a JSONL format (JSON lines), each line is an individual JSON object
        flattened_data = []
        for line in f:
            data = json.loads(line)
            reports = data.get("reports", [])
            for report in reports:
                flattened_data.append(flatten_json(report))

    # Collect all field names (headers)
    all_fieldnames = set()
    for item in flattened_data:
        all_fieldnames.update(item.keys())

    # Write the CSV file
    with open(csv_file, "w", newline="") as csvf:
        # Create a CSV DictWriter and write headers
        csv_writer = csv.DictWriter(csvf, fieldnames=sorted(all_fieldnames))
        csv_writer.writeheader()

        # Write the rows, filling missing fields with 'N/A'
        for row in flattened_data:
            csv_writer.writerow({field: row.get(field, 'N/A') for field in all_fieldnames})

    print(f"Data has been written to {csv_file}")

# Main script
json_file = "./wheat.json"  # The JSONL file you uploaded
csv_file = "./wheat.csv"   # Output CSV file

# Convert JSON to CSV
json_to_csv(json_file, csv_file)