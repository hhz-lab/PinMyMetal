import psycopg2
from psycopg2 import sql

# Connect to PostgreSQL
conn = psycopg2.connect(
    dbname="postgres",          # Use the default database
    user="myuser",               # PostgreSQL username
    password="123456",       # PostgreSQL password
    host="postgresql_container", # PostgreSQL container name as host
    port=5432                    # Default PostgreSQL port
)

# Enable autocommit to execute CREATE DATABASE
conn.autocommit = True

# Create a cursor object to execute commands
cur = conn.cursor()

# Create the new database
new_db_name = 'newdatabase'
cur.execute(sql.SQL("CREATE DATABASE {}").format(sql.Identifier(new_db_name)))

# dropdb 
cur.execute(sql.SQL("DROP DATABASE IF EXISTS {}").format(sql.Identifier(new_db_name)))


conn.commit()

# Close the cursor and connection
cur.close()
conn.close()

# Print confirmation message
print(f"Database '{new_db_name}' created successfully!")

