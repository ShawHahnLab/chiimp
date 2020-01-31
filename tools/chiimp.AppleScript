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
	set chiimp_autoclose to system attribute "CHIIMP_AUTOCLOSE"
	-- Build up a command string that incorporates the full path to the
	-- .command file, the config file argument, and the special environment variable
	-- CHIIMP_AUTOCLOSE (this isn't automatically exported to the Terminal
	-- call but it works to just build it into the string here.)
	set UnixPath to quoted form of POSIX path of ((path to me as text) & "::")
	set runCmd to UnixPath & "/exec/chiimp.command " & quoted form of POSIX path of this_item
	set runCmd to "CHIIMP_AUTOCLOSE=" & chiimp_autoclose & " " & runCmd & "; exit"
	-- Send the command string to the Terminal app, and wait for it to
	-- complete.
	tell application "Terminal"
		activate
		-- This is the cleanest way I could find to have the
		-- AppleScript wait for the shell script to finish.
		-- https://stackoverflow.com/a/17225439
		set w to do script runCmd
		repeat
			delay 0.1
			if not busy of w then exit repeat
		end repeat
	end tell
	-- Now use the CHIIMP_AUTOCLOSE environment variable here to decide if
	-- we should close the terminal window
	if (chiimp_autoclose is "yes") then
		-- https://stackoverflow.com/a/57415817
		tell application "Terminal" to close (get window 1)
	end if
end process_item
