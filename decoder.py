from cryptography.fernet import Fernet


from Crypto.PublicKey import RSA
from Crypto.Signature import PKCS1_v1_5
from base64 import b64decode



if __name__ == '__main__':
    code = 'L1UdEBQXCh1BFBluUkkCBREOGhUfGXMRAQkbEQ4JR1YedEhOQhIHGwtXXlwwVUJFUBEJCF1BTSdV Tl9XUwYAUUFcMBsMCRJTQ04VUlo8GwsTEhkKAEYUGW5SSRAZGAANWVZdc15OQgUVDQxbR0pzUlRF UAcOCFcUFXRVCAoYU09UEhROPRxPQgo= '
    key = 'trnewton239'

    rsa_key = RSA.importKey(key)
    cipher = PKCS1_v1_5.new(rsa_key)
    raw_cipher_data = b64decode(code)
    phn = cipher.decrypt(raw_cipher_data)

    print(phn)