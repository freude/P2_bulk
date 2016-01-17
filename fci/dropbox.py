import dropbox
import os
import sys
import webbrowser
 
from configobj import ConfigObj
 
########################################################################
class DropObj(object):
    """
    Dropbox object that can access your dropbox folder,
    as well as download and upload files to dropbox
    """
 
    #----------------------------------------------------------------------
    def __init__(self, filename=None, path='/'):
        """Constructor"""
        self.base_path = os.path.dirname(os.path.abspath(__file__))
        self.filename = filename
        self.path = path
        self.client = None
 
        config_path = os.path.join(self.base_path, "config.ini")
        if os.path.exists(config_path):
            try:
                cfg = ConfigObj(config_path)
            except IOError:
                print "ERROR opening config file!"
                sys.exit(1)
            self.cfg_dict = cfg.dict()
        else:
            print "ERROR: config.ini not found! Exiting!"
            sys.exit(1)
 
        self.connect()
 
    #----------------------------------------------------------------------
    def connect(self):
        """
        Connect and authenticate with dropbox
        """
        app_key = self.cfg_dict["key"]
        app_secret = self.cfg_dict["secret"]
 
        access_type = "dropbox"
        session = dropbox.session.DropboxSession(app_key,
                                                 app_secret,
                                                 access_type)
 
        request_token = session.obtain_request_token()
 
        url = session.build_authorize_url(request_token)
        msg = "Opening %s. Please make sure this application is allowed before continuing."
        print msg % url
        webbrowser.open(url)
        raw_input("Press enter to continue")
        access_token = session.obtain_access_token(request_token)
 
        self.client = dropbox.client.DropboxClient(session)
 
    #----------------------------------------------------------------------
    def download_file(self, filename=None, outDir=None):
        """
        Download either the file passed to the class or the file passed
        to the method
        """
 
        if filename:
            fname = filename
            f, metadata = self.client.get_file_and_metadata("/" + fname)
        else:
            fname = self.filename
            f, metadata = self.client.get_file_and_metadata("/" + fname)
 
        if outDir:
            dst = os.path.join(outDir, fname)
        else:
            dst = fname
 
        with open(fname, "w") as fh:
            fh.write(f.read())
 
        return dst, metadata
 
    #----------------------------------------------------------------------
    def get_account_info(self):
        """
        Returns the account information, such as user's display name,
        quota, email address, etc
        """
        return self.client.account_info()
 
    #----------------------------------------------------------------------
    def list_folder(self, folder=None):
        """
        Return a dictionary of information about a folder
        """
        if folder:
            folder_metadata = self.client.metadata(folder)
        else:
            folder_metadata = self.client.metadata("/")
        return folder_metadata
 
    #----------------------------------------------------------------------
    def upload_file(self):
        """
        Upload a file to dropbox, returns file info dict
        """
        try:
            with open(self.filename) as fh:
                path = os.path.join(self.path, self.filename)
                res = self.client.put_file(path, fh)
                print "uploaded: ", res
        except Exception, e:
            print "ERROR: ", e
 
        return res
 
if __name__ == "__main__":
    drop = DropObj("somefile.txt")
