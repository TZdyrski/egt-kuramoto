def setup_plutoserver():
  return {
    "command": ["julia", "--optimize=0", "--project=binder", "--heap-size-hint=1G", "binder/start_server.jl"],
    "environment": {
      "JSP_PORT": "{port}",
    },
    "timeout": 60,
    "launcher_entry": {
        "title": "Pluto.jl",
    },
  }
