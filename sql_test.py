import requests
import pandas as pd
import tabulate


response = submit_sql_query(142.001166,14.121277)

df,_ = response_to_dataframe(response)
print(tabulate.tabulate(df.head(), headers='keys', tablefmt='psql'))
print(_)

# df = response_to_dataframe(response)

# if df is not None:
#     print("\nDataFrame successfully created:")
#     print(df.head())
#     print(f"\nDataFrame shape: {df.shape}")
# else:
#     print("\nFailed to create DataFrame.")
