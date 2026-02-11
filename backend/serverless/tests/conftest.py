"""
Shared test configuration for backend/serverless tests.

Adds the serverless directory to sys.path so tests can import modules
with bare imports (e.g., `from handler import ...`) regardless of
which subdirectory the test lives in.
"""
import sys
import os

# Add the serverless directory (parent of tests/) to Python path
serverless_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if serverless_dir not in sys.path:
    sys.path.insert(0, serverless_dir)

# Common test output directory (gitignored)
TEST_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
os.makedirs(TEST_OUTPUT_DIR, exist_ok=True)

# Common fixtures directory
FIXTURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
