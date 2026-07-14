"""One module per workflow step (IMPLEMENTATION_BLUEPRINT.md Sec 8.11).
Each panel's build(parent, state) reads/writes gui.state.AppState and
calls only temptation.{pipeline,plotting,export,config,dataio,discovery,
compat} -- never a metric or segmentation function directly."""
