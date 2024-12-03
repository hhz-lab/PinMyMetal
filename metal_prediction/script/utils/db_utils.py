import psycopg2
import os

# Database connection parameters (can be set through environment variables)
DB_HOST = os.getenv("DATABASE_HOST", "postgresql_container")  # Default to PostgreSQL container
DB_PORT = os.getenv("DATABASE_PORT", 5432)                   # Default PostgreSQL port
DB_USER = os.getenv("DATABASE_USER", "myuser")               # Default PostgreSQL user
DB_PASSWORD = os.getenv("DATABASE_PASSWORD", "123456")       # Default PostgreSQL password

# Function to create a connection to the database
def create_connection(dbname):
    try:
        conn = psycopg2.connect(
            dbname=dbname,
            user=DB_USER,
            password=DB_PASSWORD,
            host=DB_HOST,
            port=DB_PORT
        )
        return conn
    except Exception as e:
        print(f"Error connecting to database {dbname}: {e}")
        raise

'''
# Function to execute SQL queries
def execute_query(dbname, query, params=None):
    conn = None
    try:
        conn = create_connection(dbname)
        cur = conn.cursor()
        cur.execute(query, params)
        result = cur.fetchall()
        cur.close()
        conn.commit()
        return result
    except Exception as e:
        print(f"Error executing query: {e}")
        return None
    finally:
        if conn:
            conn.close()
'''
