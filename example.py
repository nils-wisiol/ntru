from ciphers import StreamlinedNTRUPrime

# choose your parameters
p, q, w = 761, 4591, 286

print('Streamlined NTRU Prime Example for', f'p={p}, q={q}, w={w}')
print('-' * 50)
cipher = StreamlinedNTRUPrime(p, q, w, seed=1337)

print('Generating key pair ... ')
pk, sk = cipher.generate_keys()

print('En/decrypting...')
message = cipher.random_small_poly(w, None, cipher.modulus_r)
assert message == cipher.decrypt(cipher.encrypt(message, pk), sk), 'En/decryption failed.'

print('Successfully en/decrypted.')
