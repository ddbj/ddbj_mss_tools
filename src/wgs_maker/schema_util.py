import urllib.request
import json
from jsonschema import validate, ValidationError, Draft202012Validator, validators
import os
import logging
jsaon_schema_url = "https://raw.githubusercontent.com/ddbj/template_generator_api/main/src/dev_schemas/MSS_COMMON_template.json"
script_dir = os.path.dirname(os.path.abspath(__file__))
schema_local_filepath = os.path.join(script_dir, "MSS_COMMON_template.json")

logger = logging.getLogger(__name__)


def get_remote_schema():
    try:
        with urllib.request.urlopen(jsaon_schema_url) as response:
            data = response.read()
            json_data = json.loads(data.decode())
            return json_data    
    except urllib.error.HTTPError as http_err:
        print(f"HTTP error. Failed to retrieve a remote file: {http_err}")
    except Exception as err:
        print(f"Unexpected Error while retrieving a remote file: {err}")

def get_local_schema():
    if not os.path.exists(schema_local_filepath):
        download_json_file(jsaon_schema_url, schema_local_filepath)
    try:
        with open(schema_local_filepath, 'r', encoding='utf-8') as file:
            return json.load(file)
    except Exception as err:
        print(f"Error while loading a json file: {err}")
        return None

def download_json_file(url, local_filename):
    try:
        # URLからデータを取得する
        with urllib.request.urlopen(url) as response:
            data = response.read().decode('utf-8')

        # JSONデータをファイルに保存する
        with open(local_filename, 'w', encoding='utf-8') as file:
            file.write(data)

        print(f"The file has been downloaded to {local_filename}")
    except urllib.error.HTTPError as http_err:
        print(f"HTTPエラーが発生しました: {http_err}")
    except Exception as err:
        print(f"エラーが発生しました: {err}")

def load_json_file(filename):
    try:
        with open(filename, 'r', encoding='utf-8') as file:
            return json.load(file)
    except Exception as err:
        print(f"Error while loading a json file: {err}")
        return None
    

def validate_json(json_data, schema):
    try:
        validate(instance=json_data, schema=schema)
        print("JSONデータはスキーマに従っています。")
    except ValidationError as err:
        print(f"バリデーションエラー: {err.message}")
    except Exception as err:
        print(f"バリデーション中にエラーが発生しました: {err}")


def get_subschema_for_category(schema, category):
    """
    create json schema for the given category
    登録カテゴリの値に応じてjson schemaの余分な条件分岐を削除して簡潔にする
    """
    def _check_category(subschema):
        return category == subschema.get("if", {}).get("properties", {}).get("_trad_submission_category", {}).get("const", "")
    subschemas = [subschema for subschema in schema.get("allOf", []) if _check_category(subschema)]
    properties = {}
    required = []
    for subschema in subschemas:
        properties.update(subschema.get("then", {}).get("properties", {}))
        required += subschema.get("then", {}).get("required", [])
    del schema["allOf"]
    schema.setdefault("properties", {}).update(properties)
    schema.setdefault("required",[]).extend(required)


def set_default_validator(validator_class):
    # json schema のデフォルの値を設定する
    # see: https://stackoverflow.com/questions/41290777/trying-to-make-json-schema-validator-in-python-to-set-default-values
    validate_properties = validator_class.VALIDATORS["properties"]

    def set_defaults(validator, properties, instance, schema):
        for property, subschema in properties.items():
            # print(f"{property=}, {subschema=}, {instance=}")
            # if len(subschema) == 1 and "$ref" in subschema:
            # print(f"{property=}\n{subschema=}\n{instance=}\n======")
            if "$ref" in subschema:
                # if "definitions" not in schema:
                #     print("----")
                #     print(schema)
                #     exit()
                definitions = schema.get("definitions", {})
                definition_subchema = definitions[subschema["$ref"].split("/")[-1]]
                definition_subchema.update(subschema)
                subschema = definition_subchema
                del subschema["$ref"]
                # print("ref", subschema)
            if "default" in subschema:
                # print("-"*10)
                # print(f"Default value found: {subschema['default']}, {instance=}")
                instance.setdefault(property, subschema["default"])
                # print(f"Applied default, {instance=}")
                # print("-"*10)
            if subschema.get("type") == "object":
                set_defaults(validator, subschema.get("properties", {}), instance.setdefault(property, {}), schema)
            elif subschema.get("type") == "array" and property == "REFERENCE":
                for item in instance.setdefault(property, [{}]):
                    if subschema.get("items", {}).get("type") == "object":
                        set_defaults(validator, subschema.get("items", {}).get("properties", {}), item, schema)

        # for error in validate_properties(
        #     validator, properties, instance, schema,
        # ):
        #     yield error

    return validators.extend(
        validator_class, {"properties" : set_defaults},
    )

def set_default_to_json(json_data, schema):
    DefaultValidatingValidator = set_default_validator(Draft202012Validator)
    DefaultValidatingValidator(schema).validate(json_data)



if __name__ == "__main__":
    schema = get_remote_schema()
    # print(schema)
    get_subschema_for_category(schema, "MAG")
    # print(schema)
    category = "MAG"
    # json_data = {"_trad_submission_category": category, "DBLINK": {"biosample": ["SAMD99999"]}, "REFERENCE": [{"title": "title1", "ab_name": ["Tanizawa,Y.", "Fujisawa,T."]}], "SUBMITTER": {"contact": "YT", "phone": "1234567890"}}
    json_data = {"_trad_submission_category": category}
    set_default_to_json(json_data, schema)
    print(json_data)
