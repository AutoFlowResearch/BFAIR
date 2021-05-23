import Axios from 'axios';
import http from './httpService';

export const userslist = async () => {
  try {
    const { data } = await http.get(
      'https://jsonplaceholder.typicode.com/users',
    );
    return data;
  } catch (error) {
    if (error.response.status >= 400) throw new Error('error fetching users');
  }
};

export const getUserById = async (id) => {
  try {
    const { data } = await http.get(
      `https://jsonplaceholder.typicode.com/users/${id}`,
    );
    return data;
  } catch (error) {
    if (error.response.status >= 400) throw new Error('User not found');
  }
};

// export default { userslist };
