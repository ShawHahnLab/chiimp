# Mac OS AppleScript to support drag-and-drop of config file onto
# chiimp.command.  (Could alternatively use Automator.)
on open filelist
  repeat with i in filelist
    do shell script "chiimp.command " & quoted form of POSIX path of i
  end repeat
end open
