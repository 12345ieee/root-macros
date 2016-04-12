## Open in TBrowser

Open ROOT files from the graphical interface or terminal directly in a `TBrowser`,
leaving an open root shell behind.

### Install instructions
* Copy `root_open` somewhere in the `PATH` (e.g. `/usr/local/bin/`)
* Edit the `path` variable inside `root_open` to point back to the location of `open_in_TBrowser.cpp`
* Copy `TBrowser.desktop` in `~/.local/applications/`
* Now `TBrowser` will be a recognized application, just associate it to `.root` files and you're done

### Usage
Double-click on a .root file or launch `root_open $filename(s)`

