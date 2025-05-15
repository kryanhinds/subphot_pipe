# if [[ "9" == "10" ]]; then
echo $(date +%H%M)
if [ $(date +%H%M) -eq 1521 ]; then
    echo "Running tests..."
else
    echo "Not running tests..."
fi