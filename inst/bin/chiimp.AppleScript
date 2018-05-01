-- This AppleScript needs to be compiled into a Mac OS .app directory in order for drag-and-drop to work.
-- This may help get rid of the prompt that appears during each run:
-- https://discussions.apple.com/thread/2663278

property type_list : {} -- e.g.: {"PICT", "JPEG", "TIFF", "GIFf"} 
property extension_list : {"yml", "yaml", "txt"} -- e.g.: {"txt", "text", "jpg", "jpeg"}, NOT: {".txt", ".text", ".jpg", ".jpeg"}
property typeIDs_list : {} -- e.g.: {"public.jpeg", "public.tiff", "public.png"}

property post_alert : "Yes"

on run
	repeat
		display dialog "This droplet will process files dragged onto its icon." & linefeed & linefeed & "There is a user-settable preference for displaying an alert dialog when the droplet encounters a dragged-on item that is not a file of the type processed by the droplet." & return & return & "Post User Alert: " & (post_alert as text) buttons {"Set Prefs", "Done"} default button 2 with title "My File Processing Droplet"
		if the button returned of the result is "Set Prefs" then
			display dialog "Should this droplet post an alert dialog when items that are not files are dragged onto it?" & return & return & "Current Status: " & (post_alert as text) buttons {"No", "Yes"} default button post_alert
			if the button returned of the result is "Yes" then
				set post_alert to "Yes"
			else
				set post_alert to "No"
			end if
		else
			return "done"
		end if
	end repeat
end run

-- This droplet processes files dropped onto the applet 
on open these_items
	repeat with i from 1 to the count of these_items
		set this_item to item i of these_items
		set the item_info to info for this_item
		set this_name to the name of the item_info
		try
			set this_extension to the name extension of item_info
		on error
			set this_extension to ""
		end try
		try
			set this_filetype to the file type of item_info
		on error
			set this_filetype to ""
		end try
		try
			set this_typeID to the type identifier of item_info
		on error
			set this_typeID to ""
		end try
		if (folder of the item_info is false) and (alias of the item_info is false) and ((this_filetype is in the type_list) or (this_extension is in the extension_list) or (this_typeID is in typeIDs_list)) then
			-- THE ITEM IS AN IMAGE FILE AND CAN BE PROCESSED
			process_item(this_item)
		else if post_alert is "Yes" then
			display alert "PROCESSING ALERT" message "The item Ò" & this_name & "Ó is not a file that this droplet can process." buttons {"Cancel", "Continue"} default button 2 cancel button "Cancel" as warning
		end if
	end repeat
end open

on process_item(this_item)
	set UnixPath to quoted form of POSIX path of ((path to me as text) & "::")
	set runCmd to UnixPath & "/chiimp.command " & quoted form of POSIX path of this_item
	--do shell script runCmd
	tell application "Terminal"
		activate
		do shell script runCmd
	end tell
end process_item
