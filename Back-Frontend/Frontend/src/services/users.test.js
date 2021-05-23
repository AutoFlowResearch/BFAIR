import { userslist } from './users';

describe('apiCalls', () => {
  it('should return users array of object', async () => {
    const data = await userslist();
    expect(data.length).toEqual(10);
  });

  it('should return error message on fail', async () => {
    try {
      await userslist();
    } catch (e) {
      expect(e.message).toEqual('Error fetching users');
    }
  });
});
