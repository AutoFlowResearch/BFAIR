import http from './httpService';

export const postsList = async () => {
  try {
    const { data } = await http.get(
      'https://jsonplaceholder.typicode.com/posts',
    );
    return data.filter((post) => post.userId === 1);
  } catch (error) {
    if (error.response.status >= 400) throw new Error('error fetching posts');
  }
};

export const getPostById = async (id) => {
  try {
    const { data } = await http.get(
      `https://jsonplaceholder.typicode.com/posts/${id}`,
    );
    return data;
  } catch (error) {
    if (error.response.status >= 400) throw new Error('Post not found');
  }
};

export const createPost = async (newpost) => {
  try {
    const { data } = await http.post(
      'https://jsonplaceholder.typicode.com/posts',
      newpost,
    );
    // send response to client
    return data;
  } catch (error) {
    throw new Error('Post not created');
  }
};
