-- This AppleScript needs to be compiled into a Mac OS .app directory in order for drag-and-drop to work.
--
-- On the command line: osacompile -o chiimp.app chiimp.AppleScript
--
-- Interactively:
--
-- 1) Open in Script Editor
-- 2) File > Export...
--    * File Format: Application
--    * All options unchecked
--    * Export As: .../chiimp/bin/chiimp ("chiimp.app" if Hide Extension is unchecked)
--
-- When double-clicked the app should give a basic usage message.
-- When a .yml config file is dragged onto it, it should run chiimp (working dir relative to the config file).

on run
	repeat
		display dialog "CHIIMP will read configuration files dragged onto its icon." & linefeed & linefeed & "See the GUIDE.pdf file for more details on configuration." buttons {"Done"} default button 1 with title "CHIIMP"
		return "done"
	end repeat
end run

-- This droplet processes files dropped onto the applet
on open these_items
	repeat with i from 1 to the count of these_items
		set this_item to item i of these_items
		process_item(this_item)
	end repeat
end open

-- Note that the path used will only work on an installed
-- (non-inst-directory-containing) copy of the package.
on process_item(this_item)
	set UnixPath to quoted form of POSIX path of ((path to me as text) & "::")
	set runCmd to UnixPath & "/../exec/chiimp.command " & quoted form of POSIX path of this_item
	tell application "Terminal"
		activate
		do script runCmd
	end tell
end process_item
