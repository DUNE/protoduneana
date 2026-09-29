import ROOT as RT
def handle_header():
  read_header('gallery/ValidHandle.h')

def read_header(h):
  #Make the ROOT C++ jit compiler read the specified header.
  RT.gROOT.ProcessLine('#include "%s"' % h)


def provide_get_valid_handle(klass):
  """Make the ROOT C++ jit compiler instantiate the
     Event::getValidHandle member template for template
     parameter klass."""
  RT.gROOT.ProcessLine('template gallery::ValidHandle<%(name)s> gallery::Event::getValidHandle<%(name)s>(art::InputTag const&) const;' % {'name' : klass})


def provide_list(*classes):
  for c in classes:
    print('provided', c)
    provide_get_valid_handle(c)