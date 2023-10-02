# Shut down the warning for multiple copies of .gdbinit
set auto-load safe-path /

python
# Imports
import sys
sys.path.append("/home/github.com/TBKOSTER/python")
from pretty_printers import register_pretty_printers

# User instructions
register_pretty_printers()
end
