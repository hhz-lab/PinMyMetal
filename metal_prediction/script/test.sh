#!/bin/bash

DB_NAME="testdb"
DB_USER="myuser"
DB_PASSWORD="123456"
DB_HOST="postgresql_container"
DB_PORT=5432

export PGPASSWORD=$DB_PASSWORD

echo "Creating database $DB_NAME..."
psql -h $DB_HOST -U $DB_USER -d postgres -p $DB_PORT -c "CREATE DATABASE $DB_NAME;"

echo "Connecting to database $DB_NAME..."
psql -h $DB_HOST -U $DB_USER -d $DB_NAME -p $DB_PORT -c "SELECT current_database();"

echo "Database $DB_NAME is ready for further operations!"

